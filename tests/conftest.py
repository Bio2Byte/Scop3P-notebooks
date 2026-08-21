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
import common.http_lookup as http_lookup_module
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


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "real_backoff: keep the real HTTP backoff table, for tests that inspect it",
    )


@pytest.fixture(autouse=True)
def _no_retry_backoff(request, monkeypatch):
    """Remove the real pause between HTTP retries for the duration of a test.

    The backoff exists to be kind to a struggling upstream, which is worth 2.5 seconds in
    production and pure dead time in a suite -- it took the full run from 5s to 20s. The
    pause length is zeroed rather than the retry count reduced, so what is under test is
    still the real retry behaviour.
    """
    # A test that asserts on the real table has to see the real table, so it opts out.
    if request.node.get_closest_marker("real_backoff"):
        return
    monkeypatch.setattr(
        http_lookup_module,
        "BACKOFF_SECONDS",
        tuple(0 for _ in http_lookup_module.BACKOFF_SECONDS),
    )


@pytest.fixture(autouse=True)
def _no_vendored_assets(tmp_path_factory, monkeypatch):
    """Resolve browser libraries to their CDN URLs unless a test says otherwise.

    Whether assets are vendored is a property of the environment: the image has
    /opt/scop3p/vendor, a checkout does not, and CI may have either. Left unpinned, any test
    asserting on rendered view HTML passes or fails depending on where it runs -- which is
    how a rinalign view test started failing only inside the image.

    Tests that care about vendoring set SCOP3P_VENDOR_DIR themselves, and that wins over
    this because it is checked first.
    """
    import common.vendor as vendor_module

    empty = tmp_path_factory.mktemp("no-vendor")
    monkeypatch.delenv("SCOP3P_VENDOR_DIR", raising=False)
    monkeypatch.setattr(vendor_module, "DEFAULT_VENDOR_DIR", empty / "absent")
    monkeypatch.setattr(vendor_module, "_repo_vendor_dir", lambda: None)
