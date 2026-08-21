from __future__ import annotations

import contextlib
from datetime import UTC, datetime
import io
import logging
import os
from pathlib import Path
import sys
import tempfile
from typing import Any
import uuid
import warnings

from common.session_metadata import write_metadata


_CONFIGURED = False
_LOG_FILE_PATH: Path | None = None
_TRAIL_FILE_PATH: Path | None = None
_METADATA_PATH: Path | None = None
_SESSION_STARTED_AT: str | None = None
_LOG_FORMAT = "%(asctime)s %(levelname)s %(name)s event=%(event)s %(message)s"
_DEFAULT_LOG_DIR = Path("/var/log/scop3p_toolkit")


class _SafeExtraFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        if not hasattr(record, "event"):
            record.event = "-"
        return super().format(record)


def configure_logging() -> None:
    global _CONFIGURED, _LOG_FILE_PATH, _METADATA_PATH, _SESSION_STARTED_AT
    global _TRAIL_FILE_PATH
    if _CONFIGURED:
        return

    level_name = os.getenv("SCOP3P_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    session_started_at = datetime.now(UTC).isoformat()
    date_stamp = datetime.now(UTC).strftime("%Y%m%d_%H%M%S_%f")
    formatter = _SafeExtraFormatter(_LOG_FORMAT)

    log_dir = _resolve_log_dir()
    log_file_path = log_dir / f"scop3p_toolkit_log_{date_stamp}.log"
    trail_file_path = log_dir / f"scop3p_toolkit_trail_{date_stamp}.log"
    metadata_path = log_dir / "metadata.yml"

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(log_file_path, encoding="utf-8")
    file_handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.handlers = [stdout_handler, file_handler]

    # The experiment record also gets its own file, so a run can be handed over as a
    # standalone document instead of grepped out of a log interleaved with diagnostics.
    # It stays on the root logger too, so the combined log keeps the full picture.
    trail_handler = logging.FileHandler(trail_file_path, encoding="utf-8")
    trail_handler.setFormatter(formatter)
    trail_logger = logging.getLogger("scop3p.trail")
    trail_logger.addHandler(trail_handler)
    # The trail is the record: it must survive SCOP3P_LOG_LEVEL being raised to WARNING
    # for the noisy loggers, so its own level is pinned.
    trail_logger.setLevel(logging.INFO)

    write_metadata(
        metadata_path=metadata_path,
        log_file_path=log_file_path,
        log_dir=log_dir,
        session_started_at=session_started_at,
        trail_file_path=trail_file_path,
    )

    _LOG_FILE_PATH = log_file_path
    _TRAIL_FILE_PATH = trail_file_path
    _METADATA_PATH = metadata_path
    _SESSION_STARTED_AT = session_started_at
    _CONFIGURED = True

    logging.getLogger("scop3p.logging").info(
        "logging configured app=%s log_file=%s metadata=%s level=%s session_started_at=%s",
        os.getenv("SCOP3P_APP_NAME", "unknown"),
        log_file_path,
        metadata_path,
        logging.getLevelName(level),
        session_started_at,
        extra={"event": "startup"},
    )


def _resolve_log_dir() -> Path:
    requested = Path(os.getenv("SCOP3P_LOG_DIR", str(_DEFAULT_LOG_DIR)))
    try:
        requested.mkdir(parents=True, exist_ok=True)
        return requested
    except OSError:
        fallback = Path(tempfile.gettempdir()) / "scop3p_toolkit"
        fallback.mkdir(parents=True, exist_ok=True)
        return fallback


class _EventAdapter(logging.LoggerAdapter):
    def process(self, msg: str, kwargs: dict[str, Any]) -> tuple[str, dict[str, Any]]:
        extra = dict(kwargs.get("extra", {}))
        extra.setdefault("event", "-")
        kwargs["extra"] = extra
        return msg, kwargs


def get_logger(name: str) -> logging.LoggerAdapter:
    configure_logging()
    return _EventAdapter(logging.getLogger(name), {})


# ---------------------------------------------------------------------------
# The experiment trail
# ---------------------------------------------------------------------------
# A protocol run is an experiment, and the log should read as its record: which
# protocol was opened, what accession went in, which actions were taken in which
# order, and what each one produced. That is what makes a figure reproducible
# afterwards, and it is a different thing from debug output.
#
# Two properties make the difference between a record and a pile of lines:
#
#   * Ordering. Every step carries an incrementing number, so the sequence is
#     explicit rather than inferred from timestamps that can tie at millisecond
#     resolution.
#   * Attribution. Every step carries a short session id. One process serves many
#     browser sessions concurrently -- the same user with two tabs open is enough --
#     and without a discriminator two interleaved runs are indistinguishable.

#: Step verbs. A closed set, so a trail can be parsed and so the vocabulary stays
#: consistent across protocols instead of each app inventing its own wording.
TRAIL_ACTIONS = (
    "open",     # a protocol was opened
    "input",    # the user supplied a value
    "select",   # the user chose among options
    "click",    # the user triggered an action
    "result",   # that action produced something
    "blocked",  # the action could not run, and why
    "failed",   # the action raised
    "export",   # something left the app: a download or a written file
)


class ExperimentTrail:
    """An ordered, per-session record of one protocol run.

    Steps are emitted at INFO: this is the narrative a scientist or a reviewer reads,
    not diagnostics. Mechanical detail belongs at DEBUG on the app's own logger, and
    the two are separable by name (``scop3p.trail`` against ``scop3p.<app>``).
    """

    def __init__(self, app: str, session_id: str | None = None) -> None:
        self.app = app
        self.session_id = session_id or uuid.uuid4().hex[:8]
        self._step = 0
        self._logger = logging.getLogger("scop3p.trail")

    def step(self, action: str, detail: str, **fields: Any) -> int:
        """Record one step. Returns its number, which callers may quote to the user."""
        if action not in TRAIL_ACTIONS:
            # Not worth raising over -- losing a log line should never break a protocol --
            # but it must be visible, because an unknown verb breaks trail parsing.
            self._logger.warning(
                "trail action %r is not one of %s",
                action, ", ".join(TRAIL_ACTIONS),
                extra={"event": "trail_misuse"},
            )
        level = {
            "blocked": logging.WARNING,
            "failed": logging.ERROR,
        }.get(action, logging.INFO)

        self._step += 1
        suffix = "".join(f" {key}={_render(value)}" for key, value in sorted(fields.items()))
        self._logger.log(
            level,
            "app=%s session=%s step=%s action=%s detail=%s%s",
            self.app, self.session_id, self._step, action, _render(detail), suffix,
            extra={"event": "trail"},
        )
        return self._step

    # Thin named wrappers. They exist so call sites read as the thing that happened,
    # and so the verb cannot be misspelled at 150-odd sites.
    def opened(self, protocol: str, **fields: Any) -> int:
        return self.step("open", f"opened protocol {protocol}", **fields)

    def entered(self, field: str, value: Any, **fields: Any) -> int:
        # The value goes through _render so an empty accession -- a real case, the user
        # pressing Fetch on a blank field -- reads as "= -" rather than a dangling "= ".
        return self.step("input", f"{field} = {_render(value)}", **fields)

    def selected(self, field: str, value: Any, **fields: Any) -> int:
        return self.step("select", f"{field} = {_render(value)}", **fields)

    def clicked(self, label: str, **fields: Any) -> int:
        return self.step("click", label, **fields)

    def produced(self, detail: str, **fields: Any) -> int:
        return self.step("result", detail, **fields)

    def blocked(self, reason: str, **fields: Any) -> int:
        return self.step("blocked", reason, **fields)

    def failed(self, detail: str, **fields: Any) -> int:
        return self.step("failed", detail, **fields)

    def exported(self, what: str, **fields: Any) -> int:
        return self.step("export", what, **fields)

    @property
    def steps_recorded(self) -> int:
        return self._step


def new_trail(app: str | None = None, session_id: str | None = None) -> ExperimentTrail:
    """A fresh trail for one browser session.

    Call this inside the Shiny ``server`` function, never at module scope: a
    module-level trail would be shared by every session in the process and its step
    numbering would interleave two users' runs into one unreadable sequence.
    """
    configure_logging()
    return ExperimentTrail(app or os.getenv("SCOP3P_APP_NAME", "unknown"), session_id)


def _render(value: Any) -> str:
    """Quote a value only when it needs it, so lines stay readable but parseable."""
    text = "-" if value is None else str(value)
    if text == "":
        text = "-"
    return f'"{text}"' if any(character in text for character in ' \t"=') else text


# ---------------------------------------------------------------------------
# Keeping third-party chatter out of the record
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def quiet_third_party(logger: logging.LoggerAdapter | logging.Logger, event: str = "third_party"):
    """Route a library's stdout writes and warnings into the log at DEBUG.

    b2bTools reports progress with ``print`` -- 215 call sites -- and its dependencies
    warn on every prediction: scikit-learn's ``InconsistentVersionWarning`` fires once
    per unpickled model, so a testing session interleaves the protocol's own log lines
    with the same three-line warning over and over. Worse, ``warnings`` go straight to
    stderr and so never reach the log file at all.

    Nothing is discarded. Each captured line becomes a DEBUG record, which means it is
    invisible at the default INFO level and recoverable with ``SCOP3P_LOG_LEVEL=DEBUG``
    when a prediction misbehaves.

    Limitation worth knowing: this redirects Python-level ``sys.stdout``, so output
    written by native code straight to file descriptor 1 still passes through. An
    fd-level redirect would catch that too, but it would also swallow the logging
    handler's own stream, so it is deliberately not done here.
    """
    buffer = io.StringIO()
    # Bound before the with-block: if entering the redirect raises, the finally clause
    # must not itself fail with NameError and mask the original exception.
    caught: list = []
    try:
        with contextlib.redirect_stdout(buffer), warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            caught = recorded
            yield
    finally:
        for line in buffer.getvalue().splitlines():
            if line.strip():
                logger.debug("%s stdout: %s", event, line.rstrip(), extra={"event": event})
        for warning in caught:
            logger.debug(
                "%s warning: %s: %s",
                event,
                getattr(warning.category, "__name__", warning.category),
                str(warning.message).strip().replace("\n", " "),
                extra={"event": event},
            )


def get_log_file_path() -> Path | None:
    return _LOG_FILE_PATH


def get_trail_file_path() -> Path | None:
    """The standalone experiment record for this session."""
    return _TRAIL_FILE_PATH


def get_metadata_path() -> Path | None:
    return _METADATA_PATH


def get_session_started_at() -> str | None:
    return _SESSION_STARTED_AT


def _reset_logging_for_tests() -> None:
    global _CONFIGURED, _LOG_FILE_PATH, _METADATA_PATH, _SESSION_STARTED_AT
    global _TRAIL_FILE_PATH
    for handler in logging.getLogger().handlers:
        handler.close()
    logging.getLogger().handlers = []
    trail_logger = logging.getLogger("scop3p.trail")
    for handler in trail_logger.handlers:
        handler.close()
    trail_logger.handlers = []
    trail_logger.setLevel(logging.NOTSET)
    _CONFIGURED = False
    _TRAIL_FILE_PATH = None
    _LOG_FILE_PATH = None
    _METADATA_PATH = None
    _SESSION_STARTED_AT = None
