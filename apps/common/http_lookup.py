"""One HTTP policy for the annotation lookups every protocol makes.

These calls run inside synchronous Shiny reactive effects, so an unbounded request does
not just delay the user who made it -- it blocks the worker, and every other connected
session with it. The policy is therefore deliberately impatient, and it is shared so a
protocol cannot quietly opt out of it.

Both failure modes here were observed against real upstreams while testing, not imagined:

* **A flaky connect or handshake.** ``rest.uniprot.org`` has timed out on connect, and
  dropped a TLS handshake with ``UNEXPECTED_EOF_WHILE_READING``, in both cases while a
  direct request moments earlier answered in well under a second. One retry recovers the
  run rather than leaving the user with empty dropdowns.
* **A host that accepts the connection and then stalls.** The EBI Proteins API has
  connected in 0.06s and then never answered. Only the read timeout bounds that, which is
  why it is far below what a file download would want.

File downloads are deliberately *not* routed through here: a large structure legitimately
takes time, and giving it this read timeout would break it.
"""

from __future__ import annotations

import logging
import time
from typing import Callable

import requests

#: Seconds to wait for a TCP connect (and TLS handshake). Short on purpose: an
#: unreachable host should fail fast rather than hold the UI.
CONNECT_TIMEOUT = 5

#: Seconds to wait for the response body of a metadata lookup.
READ_TIMEOUT = 20

#: Total attempts. Three, not two: with a pause between them a second retry is cheap and
#: covers a fault that lasts a second or so, while a real outage still surfaces in a few
#: seconds rather than being retried into a long stall.
ATTEMPTS = 3

#: Seconds to wait *before* each retry. Retrying instantly is the least useful thing to do
#: to a server that is already struggling or shedding load -- it doubles the instantaneous
#: burst at the exact moment the far end is asking for less. The pause also gives a
#: transient fault time to pass, which an immediate retry does not. One entry per retry, so
#: len(BACKOFF_SECONDS) + 1 == ATTEMPTS.
BACKOFF_SECONDS = (0.5, 2.0)

LOOKUP_TIMEOUT = (CONNECT_TIMEOUT, READ_TIMEOUT)


def lookup(
    url: str,
    *,
    logger: logging.Logger | logging.LoggerAdapter | None = None,
    attempts: int = ATTEMPTS,
    timeout: tuple[int, int] = LOOKUP_TIMEOUT,
    validate: Callable[[requests.Response], object] | None = None,
    sleep: Callable[[float], None] = time.sleep,
    **kwargs,
) -> requests.Response:
    """GET a metadata lookup with bounded timeouts, backoff and retries.

    ``validate`` is called on each response *inside* the retry loop, and anything it raises
    is treated as a failed attempt. That is how a truncated body gets retried: a partial
    read is a transport fault, but it only becomes visible when something tries to parse
    the body, which happens after the request has already "succeeded". Observed for real
    against UniProt -- ``Expecting ',' delimiter: line 1 column 88423`` -- and without this
    it escaped the retry entirely and surfaced as a failed lookup.

    Raises the last error if every attempt fails: retrying must never turn a real outage
    into a silently empty result, because an empty dropdown reads to the user as "this
    protein has no structures".
    """
    total = max(1, attempts)
    last_error: Exception | None = None
    for attempt in range(1, total + 1):
        try:
            response = requests.get(url, timeout=timeout, **kwargs)
            if validate is not None:
                validate(response)
            return response
        except (requests.RequestException, ValueError) as error:
            # ValueError covers json.JSONDecodeError, which subclasses it -- that is the
            # truncated-body case above, not a programming error being swallowed.
            last_error = error
            if logger is not None:
                logger.warning(
                    "lookup failed attempt=%s/%s url=%s error=%s",
                    attempt,
                    total,
                    url,
                    error,
                    extra={"event": "http_retry"},
                )
            if attempt < total:
                # No pause after the final attempt: there is nothing left to wait for.
                index = min(attempt - 1, len(BACKOFF_SECONDS) - 1)
                sleep(BACKOFF_SECONDS[index])
    raise last_error  # type: ignore[misc]


def json_body_validator(response: requests.Response) -> None:
    """Validate that a response which *claims* to be JSON carries a complete body.

    Pass as ``validate`` where the caller needs the ``Response`` itself -- to check a status
    code or a content type -- but still wants a truncated body retried.

    Two conditions, and both matter:

    * **2xx only.** A 404 or a 500 is an *answer*. Retrying it wastes the user's time and
      adds load to a server that already replied, and the caller's status handling is what
      should deal with it.
    * **A JSON content type only.** Scop3P's single-page-app catch-all answers 200 with
      *HTML* when an endpoint has moved. That is a permanent condition, not a truncated
      read: retrying it is pointless, and parsing it here would replace the client's
      actionable "the endpoint has most likely moved" message with a bare JSON error.

    What is left is exactly the transient case: a 200, declared JSON, body cut short.

    The body is parsed once here and again by the caller. That is a millisecond on payloads
    this size, and the alternative -- threading the parsed value back out -- buys nothing.
    """
    if not 200 <= response.status_code < 300:
        return
    content_type = (response.headers.get("content-type") or "").split(";")[0].strip()
    if content_type.lower() != "application/json":
        return
    response.json()


def lookup_json(
    url: str,
    *,
    logger: logging.Logger | logging.LoggerAdapter | None = None,
    **kwargs,
):
    """GET and parse JSON, with the parse inside the retry loop.

    Use this rather than ``lookup(...).json()`` wherever only the payload is wanted: the
    latter parses *after* the retry loop has finished, so a truncated body cannot be
    retried. The payload is parsed exactly once.
    """
    parsed: dict[str, object] = {}

    def validate(response: requests.Response) -> None:
        parsed["payload"] = response.json()

    lookup(url, logger=logger, validate=validate, **kwargs)
    return parsed["payload"]
