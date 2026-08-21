"""Shared fixtures.

The cache is process-wide by design, which is exactly what makes it useful in production
and a hazard in a test suite: without clearing it between tests, one test's cached
response satisfies another's "and now it fetches" assertion, and the suite becomes
order-dependent. That failure is intermittent and extremely annoying to trace, so it is
prevented here rather than remembered at each call site.
"""

from __future__ import annotations

import pytest

import common.cache as cache_module
from common.cache import REGISTRY, clear_all


@pytest.fixture(autouse=True)
def _isolate_lookup_caches():
    """Every test starts with empty caches and clean counters."""
    clear_all()
    for cache in REGISTRY.values():
        cache.stats.__init__()  # reset counters without rebinding the object
    yield
    clear_all()


@pytest.fixture(autouse=True)
def _isolate_structure_cache(tmp_path_factory, monkeypatch):
    """Give each test its own structure-file directory.

    The shared directory is process-wide and lives under the system temp dir, so it
    outlives a test run. A test asserting "and then it downloads" would silently pass on
    a file left behind by an earlier run -- the same order-dependence the memo cache
    needed protecting from, but persisting across runs as well.
    """
    directory = tmp_path_factory.mktemp("structures")
    monkeypatch.setenv("SCOP3P_STRUCTURE_CACHE_DIR", str(directory))
    monkeypatch.setattr(cache_module, "_STRUCTURE_DIR", None, raising=False)
    yield directory
    monkeypatch.setattr(cache_module, "_STRUCTURE_DIR", None, raising=False)
