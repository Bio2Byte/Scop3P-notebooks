from __future__ import annotations

import logging
import os
import sys
from typing import Any


_CONFIGURED = False


def configure_logging() -> None:
    global _CONFIGURED
    if _CONFIGURED:
        return

    level_name = os.getenv("SCOP3P_LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)

    logging.basicConfig(
        level=level,
        stream=sys.stdout,
        format="%(asctime)s %(levelname)s %(name)s event=%(event)s %(message)s",
    )
    _CONFIGURED = True


class _EventAdapter(logging.LoggerAdapter):
    def process(self, msg: str, kwargs: dict[str, Any]) -> tuple[str, dict[str, Any]]:
        extra = dict(kwargs.get("extra", {}))
        extra.setdefault("event", "-")
        kwargs["extra"] = extra
        return msg, kwargs


def get_logger(name: str) -> logging.LoggerAdapter:
    configure_logging()
    return _EventAdapter(logging.getLogger(name), {})

