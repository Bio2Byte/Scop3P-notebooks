from __future__ import annotations

from datetime import UTC, datetime
import logging
import os
from pathlib import Path
import sys
import tempfile
from typing import Any

from common.session_metadata import write_metadata


_CONFIGURED = False
_LOG_FILE_PATH: Path | None = None
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
    if _CONFIGURED:
        return

    level_name = os.getenv("SCOP3P_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)
    session_started_at = datetime.now(UTC).isoformat()
    date_stamp = datetime.now(UTC).strftime("%Y%m%d_%H%M%S_%f")
    formatter = _SafeExtraFormatter(_LOG_FORMAT)

    log_dir = _resolve_log_dir()
    log_file_path = log_dir / f"scop3p_toolkit_log_{date_stamp}.log"
    metadata_path = log_dir / "metadata.yml"

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(log_file_path, encoding="utf-8")
    file_handler.setFormatter(formatter)

    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.handlers = [stdout_handler, file_handler]

    write_metadata(
        metadata_path=metadata_path,
        log_file_path=log_file_path,
        log_dir=log_dir,
        session_started_at=session_started_at,
    )

    _LOG_FILE_PATH = log_file_path
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


def log_action_button_click(logger: logging.LoggerAdapter, button_id: str, click_count: Any) -> None:
    logger.info(
        "action_button clicked button=%s click_count=%s",
        button_id,
        click_count,
        extra={"event": "action_button_click"},
    )


def get_log_file_path() -> Path | None:
    return _LOG_FILE_PATH


def get_metadata_path() -> Path | None:
    return _METADATA_PATH


def get_session_started_at() -> str | None:
    return _SESSION_STARTED_AT


def _reset_logging_for_tests() -> None:
    global _CONFIGURED, _LOG_FILE_PATH, _METADATA_PATH, _SESSION_STARTED_AT
    for handler in logging.getLogger().handlers:
        handler.close()
    logging.getLogger().handlers = []
    _CONFIGURED = False
    _LOG_FILE_PATH = None
    _METADATA_PATH = None
    _SESSION_STARTED_AT = None
