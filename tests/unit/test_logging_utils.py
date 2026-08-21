from __future__ import annotations

import contextlib
import io
import logging
import re

from common.logging_utils import (
    _SafeExtraFormatter,
    _EventAdapter,
    _reset_logging_for_tests,
    configure_logging,
    get_log_file_path,
    get_metadata_path,
    new_trail,
)


@contextlib.contextmanager
def _capture(logger: logging.Logger):
    """Collect records emitted by one logger, without touching global config."""
    records: list[logging.LogRecord] = []

    class _Collector(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record)

    handler = _Collector()
    previous_handlers, previous_level = logger.handlers, logger.level
    logger.handlers = [handler]
    logger.setLevel(logging.DEBUG)
    try:
        yield records
    finally:
        logger.handlers = previous_handlers
        logger.setLevel(previous_level)


def test_safe_extra_formatter_defaults_missing_event() -> None:
    stream = io.StringIO()
    handler = logging.StreamHandler(stream)
    handler.setFormatter(_SafeExtraFormatter("%(levelname)s event=%(event)s %(message)s"))

    logger = logging.getLogger("test.logging_utils")
    logger.handlers = [handler]
    logger.propagate = False
    logger.setLevel(logging.INFO)

    logger.info("hello")

    assert "INFO event=- hello" in stream.getvalue()


def test_a_click_is_recorded_as_an_ordered_trail_step() -> None:
    """Replaces the old log_action_button_click helper.

    A raw click count told you a button was pressed; the trail says which action it was,
    where it sits in the sequence, and which session it belongs to.
    """
    logger = logging.getLogger("scop3p.trail")
    with _capture(logger) as records:
        trail = new_trail("mutation-effect", session_id="sess0001")
        trail.clicked("Fetch + Predict WT")

    message = records[-1].getMessage()
    assert "action=click" in message
    assert "Fetch + Predict WT" in message
    assert "step=1" in message
    assert "session=sess0001" in message



def test_configure_logging_mirrors_records_to_log_file(tmp_path, monkeypatch) -> None:
    _reset_logging_for_tests()
    monkeypatch.setenv("SCOP3P_LOG_DIR", str(tmp_path))
    monkeypatch.setenv("SCOP3P_LOG_LEVEL", "INFO")
    monkeypatch.setenv("SCOP3P_APP_NAME", "test-app")

    try:
        configure_logging()
        logger = logging.getLogger("test.file_logging")
        logger.info("hello file", extra={"event": "unit"})

        log_file = get_log_file_path()
        assert log_file is not None
        assert log_file.parent == tmp_path
        assert re.fullmatch(r"scop3p_toolkit_log_\d{8}_\d{6}_\d{6}\.log", log_file.name)

        contents = log_file.read_text(encoding="utf-8")
        assert "INFO test.file_logging event=unit hello file" in contents
        assert "INFO scop3p.logging event=startup logging configured app=test-app" in contents
    finally:
        _reset_logging_for_tests()


def test_configure_logging_is_idempotent(tmp_path, monkeypatch) -> None:
    _reset_logging_for_tests()
    monkeypatch.setenv("SCOP3P_LOG_DIR", str(tmp_path))

    try:
        configure_logging()
        first_log_file = get_log_file_path()
        first_handlers = list(logging.getLogger().handlers)

        configure_logging()

        assert get_log_file_path() == first_log_file
        assert logging.getLogger().handlers == first_handlers
        assert len(logging.getLogger().handlers) == 2
    finally:
        _reset_logging_for_tests()


def test_configure_logging_writes_metadata(tmp_path, monkeypatch) -> None:
    _reset_logging_for_tests()
    monkeypatch.setenv("SCOP3P_LOG_DIR", str(tmp_path))
    monkeypatch.setenv("SCOP3P_APP_NAME", "metadata-test")
    monkeypatch.setenv("SCOP3P_IMAGE_VERSION", "v1.2.3")
    monkeypatch.setenv("SCOP3P_IMAGE_REVISION", "abc123")
    monkeypatch.setenv("SCOP3P_IMAGE_CREATED", "2026-05-19T00:00:00Z")

    try:
        configure_logging()

        metadata_path = get_metadata_path()
        log_file = get_log_file_path()
        assert metadata_path == tmp_path / "metadata.yml"
        assert log_file is not None

        contents = metadata_path.read_text(encoding="utf-8")
        assert 'name: "metadata-test"' in contents
        assert 'version: "v1.2.3"' in contents
        assert 'revision: "abc123"' in contents
        assert f'log_directory: "{tmp_path}"' in contents
        assert f'log_file: "{log_file}"' in contents
        assert "dependencies:" in contents
        assert "external_tools:" in contents
    finally:
        _reset_logging_for_tests()
