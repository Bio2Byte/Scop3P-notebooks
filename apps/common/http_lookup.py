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

import requests

#: Seconds to wait for a TCP connect (and TLS handshake). Short on purpose: an
#: unreachable host should fail fast rather than hold the UI.
CONNECT_TIMEOUT = 5

#: Seconds to wait for the response body of a metadata lookup.
READ_TIMEOUT = 20

#: Total attempts, so exactly one retry. Enough for an intermittent fault, while a real
#: outage still surfaces promptly instead of being retried into a long stall.
ATTEMPTS = 2

LOOKUP_TIMEOUT = (CONNECT_TIMEOUT, READ_TIMEOUT)


def lookup(
    url: str,
    *,
    logger: logging.Logger | logging.LoggerAdapter | None = None,
    attempts: int = ATTEMPTS,
    timeout: tuple[int, int] = LOOKUP_TIMEOUT,
    **kwargs,
) -> requests.Response:
    """GET a metadata lookup with bounded timeouts and one retry.

    Raises the last transport error if every attempt fails: retrying must never turn a
    real outage into a silently empty result, because an empty dropdown reads to the user
    as "this protein has no structures".
    """
    last_error: requests.RequestException | None = None
    for attempt in range(1, max(1, attempts) + 1):
        try:
            return requests.get(url, timeout=timeout, **kwargs)
        except requests.RequestException as error:
            last_error = error
            if logger is not None:
                logger.warning(
                    "lookup failed attempt=%s/%s url=%s error=%s",
                    attempt,
                    attempts,
                    url,
                    error,
                    extra={"event": "http_retry"},
                )
    raise last_error  # type: ignore[misc]
