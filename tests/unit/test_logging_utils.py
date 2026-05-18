from __future__ import annotations

import io
import logging

from common.logging_utils import _SafeExtraFormatter


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
