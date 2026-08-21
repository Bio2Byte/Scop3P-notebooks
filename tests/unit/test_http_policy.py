"""The shared HTTP policy: bounded timeouts, backoff, and what counts as retryable.

Both behaviours here were added in response to observed failures against real upstreams,
not anticipated ones:

* ``rest.uniprot.org`` returning a **truncated JSON body** -- "Expecting ',' delimiter:
  line 1 column 88423". A partial read is a transport fault, but it only becomes visible
  when something parses the body, which is after the request has already "succeeded", so it
  escaped the retry entirely.
* The same host dropping TLS handshakes and closing connections mid-response, repeatedly.
  Retrying those instantly doubles the burst on a server that is evidently struggling.
"""

from __future__ import annotations

import json

import pytest
import requests

from common import http_lookup as policy
from common.http_lookup import json_body_validator, lookup, lookup_json


class _Response:
    def __init__(self, payload=None, *, status_code=200, text=None, content_type="application/json"):
        self._payload = payload
        self.status_code = status_code
        self.text = text if text is not None else json.dumps(payload)
        self.headers = {"content-type": content_type}

    def json(self):
        return json.loads(self.text)

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")


@pytest.fixture
def sleeps():
    """Record the pauses a lookup asks for, without taking them."""
    return []


# --------------------------------------------------------------------------------------
# Backoff
# --------------------------------------------------------------------------------------


def test_the_pause_grows_between_attempts(monkeypatch, sleeps) -> None:
    """Retrying instantly is the least useful thing to do to a struggling server."""
    monkeypatch.setattr(policy, "BACKOFF_SECONDS", (0.5, 2.0))
    monkeypatch.setattr(
        "common.http_lookup.requests.get",
        lambda *a, **k: (_ for _ in ()).throw(requests.ConnectTimeout("nope")),
    )
    with pytest.raises(requests.ConnectTimeout):
        lookup("https://example.invalid", sleep=sleeps.append)
    assert sleeps == [0.5, 2.0]


def test_there_is_no_pause_after_the_final_attempt(monkeypatch, sleeps) -> None:
    """Nothing is waiting on it; the delay would be pure dead time for the user."""
    monkeypatch.setattr(policy, "BACKOFF_SECONDS", (0.5, 2.0))
    monkeypatch.setattr(
        "common.http_lookup.requests.get",
        lambda *a, **k: (_ for _ in ()).throw(requests.ConnectTimeout("nope")),
    )
    with pytest.raises(requests.ConnectTimeout):
        lookup("https://example.invalid", sleep=sleeps.append)
    assert len(sleeps) == policy.ATTEMPTS - 1


def test_a_success_never_pauses(monkeypatch, sleeps) -> None:
    monkeypatch.setattr(
        "common.http_lookup.requests.get", lambda *a, **k: _Response({"ok": True})
    )
    lookup("https://example.test", sleep=sleeps.append)
    assert sleeps == []


def test_one_pause_when_the_second_attempt_succeeds(monkeypatch, sleeps) -> None:
    attempts = {"n": 0}

    def flaky(*args, **kwargs):
        attempts["n"] += 1
        if attempts["n"] == 1:
            raise requests.ConnectionError("closed without response")
        return _Response({"ok": True})

    monkeypatch.setattr("common.http_lookup.requests.get", flaky)
    lookup("https://example.test", sleep=sleeps.append)
    assert len(sleeps) == 1


@pytest.mark.real_backoff
def test_the_backoff_table_matches_the_attempt_count() -> None:
    """One pause per retry; a mismatch would silently drop or repeat a delay."""
    assert len(policy.BACKOFF_SECONDS) == policy.ATTEMPTS - 1
    assert list(policy.BACKOFF_SECONDS) == sorted(policy.BACKOFF_SECONDS)
    assert policy.BACKOFF_SECONDS[0] > 0


@pytest.mark.real_backoff
def test_a_real_outage_still_surfaces_quickly() -> None:
    """The point is politeness, not patience: the user is waiting."""
    worst_case = sum(policy.BACKOFF_SECONDS) + policy.ATTEMPTS * policy.CONNECT_TIMEOUT
    assert worst_case < 30, f"a dead host would hold the user for {worst_case}s"


# --------------------------------------------------------------------------------------
# Truncated bodies
# --------------------------------------------------------------------------------------

TRUNCATED = '{"uniProtKBCrossReferences": [{"database": "PDB", "id": "2IV'


def test_a_truncated_json_body_is_retried(monkeypatch, sleeps) -> None:
    """The failure that escaped the retry: it only shows up at parse time."""
    attempts = {"n": 0}

    def flaky(*args, **kwargs):
        attempts["n"] += 1
        if attempts["n"] == 1:
            return _Response(text=TRUNCATED)
        return _Response({"complete": True})

    monkeypatch.setattr("common.http_lookup.requests.get", flaky)
    payload = lookup_json("https://example.test", sleep=sleeps.append)
    assert payload == {"complete": True}
    assert attempts["n"] == 2


def test_a_persistently_truncated_body_raises(monkeypatch, sleeps) -> None:
    """Retrying must not turn a broken upstream into a silently empty result."""
    monkeypatch.setattr(
        "common.http_lookup.requests.get", lambda *a, **k: _Response(text=TRUNCATED)
    )
    with pytest.raises(ValueError):
        lookup_json("https://example.test", sleep=sleeps.append)


def test_the_payload_is_parsed_only_once(monkeypatch, sleeps) -> None:
    """lookup_json exists so the parse happens inside the loop, not twice outside it."""
    parses = {"n": 0}

    class _Counting(_Response):
        def json(self):
            parses["n"] += 1
            return super().json()

    monkeypatch.setattr(
        "common.http_lookup.requests.get", lambda *a, **k: _Counting({"ok": True})
    )
    assert lookup_json("https://example.test", sleep=sleeps.append) == {"ok": True}
    assert parses["n"] == 1


# --------------------------------------------------------------------------------------
# What must NOT be retried
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("status", [400, 404, 410, 500, 503])
def test_an_http_error_is_not_retried(monkeypatch, sleeps, status) -> None:
    """A status code is an answer. Retrying it wastes the user's time and adds load."""
    attempts = {"n": 0}

    def responder(*args, **kwargs):
        attempts["n"] += 1
        return _Response({"error": True}, status_code=status)

    monkeypatch.setattr("common.http_lookup.requests.get", responder)
    lookup("https://example.test", validate=json_body_validator, sleep=sleeps.append)
    assert attempts["n"] == 1
    assert sleeps == []


def test_an_html_body_on_200_is_not_retried(monkeypatch, sleeps) -> None:
    """Scop3P's single-page-app catch-all answers 200 with HTML when a path has moved.

    That is permanent, not a truncated read. Retrying it is pointless, and parsing it here
    would replace the client's actionable "the endpoint has most likely moved" message with
    a bare JSON error -- which is exactly what happened before the content-type condition.
    """
    attempts = {"n": 0}

    def responder(*args, **kwargs):
        attempts["n"] += 1
        return _Response(text="<!doctype html><html></html>", content_type="text/html")

    monkeypatch.setattr("common.http_lookup.requests.get", responder)
    lookup("https://example.test", validate=json_body_validator, sleep=sleeps.append)
    assert attempts["n"] == 1, "the SPA catch-all was retried"
    assert sleeps == []


def test_the_validator_ignores_a_non_json_content_type() -> None:
    """FASTA and PDB text must pass straight through."""
    json_body_validator(_Response(text=">sp|P1|\nACDE", content_type="text/plain"))


def test_the_validator_rejects_a_truncated_json_body() -> None:
    with pytest.raises(ValueError):
        json_body_validator(_Response(text=TRUNCATED))


def test_the_validator_accepts_a_complete_json_body() -> None:
    json_body_validator(_Response({"complete": True}))
