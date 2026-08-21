"""A process-wide memo for external lookups, shared by every protocol.

Six protocols run in one process behind the portal, and they ask the same upstreams for
the same things: the UniProt sequence for an accession is fetched by Mutation Effect and
again by Structure Visualisation, Scop3P modifications by four of them. Each session
repeated all of it. Against upstreams that fail intermittently -- UniProt has dropped TLS
handshakes and the EBI Proteins API has stalled mid-response -- every avoidable request is
another chance to fail in front of the user.

Design notes, in rough order of how much trouble each would cause if got wrong:

**Failures are never cached.** This is the one that matters. The upstreams here fail
transiently, so caching an exception would turn a one-off blip into a dead protocol: the
Fetch button would keep "failing" from cache with no request ever made, and the only cure
would be a restart. A raised exception leaves no entry behind, so a retry is a real retry.

**The key ignores ``self``.** These are service methods and each browser session builds its
own service instance, so keying on the bound instance would give every session a private
cache and share nothing -- the entire point missed, invisibly, while still looking like a
working cache.

**Mutable results are copied on the way out.** A DataFrame handed to two sessions is one
object; if either mutates it the other silently sees the change. Copying on return costs a
little memory and removes a class of bug that would be very hard to trace back to here.

**Why not memcached.** It was worth considering and it is the wrong shape for this. The
data is a few MB per accession, it does not need to outlive the process, and there is one
process per container -- so a cache server would add a service to every image, a network
hop, serialisation, and a new failure mode, in exchange for nothing this needs. A bounded
dict in the same address space is faster and cannot be down. If the toolkit is ever scaled
to several worker processes or replicas, that is the point to revisit it: then the sharing
is genuinely cross-process and the trade changes.
"""

from __future__ import annotations

import functools
import inspect
import os
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from common.logging_utils import get_logger

LOGGER = get_logger("scop3p.common.cache")

#: How long an entry stays fresh. Sequences, PTM annotations and structure
#: cross-references change on release timescales, not within a session, so this is about
#: bounding staleness for a long-running container rather than correctness.
DEFAULT_TTL_SECONDS = 60 * 60

#: Entries per cache. Sized for JSON and small frames; structure *files* belong on disk,
#: not here. A few hundred accessions is well beyond a working session.
DEFAULT_MAX_ENTRIES = 256

#: Bio2Byte predictions are far larger per entry than a JSON lookup -- roughly 150 KB for
#: a 1100-residue protein, and proportionally more for a long one -- so they get their own
#: smaller bound. They are also the most expensive thing the toolkit computes (about 16
#: seconds for 1100 residues), which is what makes caching them worth the memory.
PREDICTION_MAX_ENTRIES = 32

#: Predictions are deterministic for a given sequence, so there is nothing to go stale
#: within the life of a process: a b2bTools upgrade arrives as a new image, which restarts
#: it. The long TTL only bounds memory in a very long-lived container.
PREDICTION_TTL_SECONDS = 24 * 60 * 60


@dataclass
class CacheStats:
    hits: int = 0
    misses: int = 0
    evictions: int = 0
    expirations: int = 0
    errors_not_cached: int = 0

    def as_dict(self) -> dict[str, int]:
        return {
            "hits": self.hits,
            "misses": self.misses,
            "evictions": self.evictions,
            "expirations": self.expirations,
            "errors_not_cached": self.errors_not_cached,
        }


@dataclass
class _Entry:
    value: Any
    stored_at: float


@dataclass
class LookupCache:
    """A bounded, expiring, thread-safe memo.

    Thread safety is required, not defensive: blocking handlers run under
    ``asyncio.to_thread``, so two sessions can be inside the same lookup at once.
    """

    name: str
    ttl_seconds: float = DEFAULT_TTL_SECONDS
    max_entries: int = DEFAULT_MAX_ENTRIES
    stats: CacheStats = field(default_factory=CacheStats)
    _entries: dict[Any, _Entry] = field(default_factory=dict)
    _lock: threading.RLock = field(default_factory=threading.RLock)

    def get(self, key: Any) -> tuple[bool, Any]:
        """``(found, value)``. A tuple rather than a sentinel so ``None`` can be cached."""
        with self._lock:
            entry = self._entries.get(key)
            if entry is None:
                self.stats.misses += 1
                return False, None
            if self._expired(entry):
                del self._entries[key]
                self.stats.expirations += 1
                self.stats.misses += 1
                return False, None
            # Refresh recency for the LRU bound.
            self._entries[key] = self._entries.pop(key)
            self.stats.hits += 1
            return True, entry.value

    def set(self, key: Any, value: Any) -> None:
        with self._lock:
            if key in self._entries:
                del self._entries[key]
            self._entries[key] = _Entry(value=value, stored_at=time.monotonic())
            while len(self._entries) > self.max_entries:
                self._entries.pop(next(iter(self._entries)))
                self.stats.evictions += 1

    def invalidate(self, predicate: Callable[[Any], bool] | None = None) -> int:
        """Drop everything, or every key the predicate accepts. Returns the count."""
        with self._lock:
            keys = [
                key for key in self._entries if predicate is None or predicate(key)
            ]
            for key in keys:
                del self._entries[key]
            return len(keys)

    def _expired(self, entry: _Entry) -> bool:
        return (
            self.ttl_seconds is not None
            and (time.monotonic() - entry.stored_at) > self.ttl_seconds
        )

    def __len__(self) -> int:
        with self._lock:
            return len(self._entries)


#: Every cache built by @memoize, so a status page or a test can inspect them.
REGISTRY: dict[str, LookupCache] = {}


def _copy_out(value: Any) -> Any:
    """Hand back a copy of anything a caller could mutate.

    Frozen dataclasses inside a list are safe individually, but the list is not: a caller
    doing ``refs.sort()`` would reorder the cached value for every other session.
    """
    if isinstance(value, pd.DataFrame):
        return value.copy()
    if isinstance(value, pd.Series):
        return value.copy()
    if isinstance(value, list):
        return list(value)
    if isinstance(value, dict):
        return dict(value)
    if isinstance(value, set):
        return set(value)
    return value


def memoize(
    *,
    name: str | None = None,
    ttl_seconds: float = DEFAULT_TTL_SECONDS,
    max_entries: int = DEFAULT_MAX_ENTRIES,
    skip_self: bool = True,
):
    """Memoize a lookup process-wide, across sessions and protocols.

    For a **method**, the bound instance is dropped from the key. That is what makes a
    service method shareable: every session constructs its own service, and keying on the
    instance would silently give each one a private cache.

    Whether to drop it is decided from the signature -- the first parameter being named
    ``self`` or ``cls`` -- not from the call. Dropping the first positional argument
    unconditionally looks equivalent and is not: applied to a plain function it discards
    the real argument, so ``fetch("P1")`` and ``fetch("P2")`` collapse onto one entry and
    the second accession quietly returns the first one's data.

    An exception is propagated and nothing is stored -- see the module docstring.
    """

    def decorate(function: Callable) -> Callable:
        cache_name = name or f"{function.__module__}.{function.__qualname__}"
        try:
            first_parameter = next(iter(inspect.signature(function).parameters), None)
        except (TypeError, ValueError):  # pragma: no cover - builtins have no signature
            first_parameter = None
        drop_receiver = skip_self and first_parameter in {"self", "cls"}
        # The name *is* the cache identity, and reusing one is how two protocols share a
        # result: Mutation Effect and Structure Visualisation both fetch the same UniProt
        # FASTA, so both declare "uniprot.sequence.fasta" and hit one store. Creating a
        # fresh cache per decorated function instead would give each its own copy while
        # still looking shared -- and would orphan all but the last from REGISTRY, so
        # clear_all() would silently miss them.
        cache = REGISTRY.get(cache_name)
        if cache is None:
            cache = LookupCache(
                name=cache_name, ttl_seconds=ttl_seconds, max_entries=max_entries
            )
            REGISTRY[cache_name] = cache

        @functools.wraps(function)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            key_args = args[1:] if (drop_receiver and args) else args
            try:
                key = (key_args, tuple(sorted(kwargs.items())))
            except TypeError:
                # An unhashable argument: run uncached rather than fail. Nothing here
                # takes one today, so this is a guard, not a routine path.
                LOGGER.debug(
                    "cache bypassed for %s: unhashable arguments",
                    cache_name,
                    extra={"event": "cache"},
                )
                return function(*args, **kwargs)

            found, value = cache.get(key)
            if found:
                LOGGER.debug(
                    "cache hit %s key=%s", cache_name, key_args, extra={"event": "cache"}
                )
                return _copy_out(value)

            try:
                result = function(*args, **kwargs)
            except Exception:
                # Deliberately not cached. A transient upstream failure must not become a
                # permanent one, or the Fetch button would stop being a retry.
                cache.stats.errors_not_cached += 1
                raise

            cache.set(key, result)
            LOGGER.debug(
                "cache store %s key=%s entries=%s",
                cache_name,
                key_args,
                len(cache),
                extra={"event": "cache"},
            )
            return _copy_out(result)

        wrapper.cache = cache  # type: ignore[attr-defined]
        return wrapper

    return decorate


def cache_report() -> dict[str, dict[str, int]]:
    """Hit/miss counters per cache, for logging or a status readout."""
    return {name: cache.stats.as_dict() for name, cache in sorted(REGISTRY.items())}


def clear_all() -> None:
    """Drop every cached entry. Used by tests, and available to an operator."""
    for cache in REGISTRY.values():
        cache.invalidate()


# ---------------------------------------------------------------------------
# Structure files
# ---------------------------------------------------------------------------
# Downloaded structures are immutable upstream artefacts: 2IVT is 2IVT for every session
# and every protocol. They were being written into each session's own working directory,
# so opening two tabs downloaded the same megabytes twice, and each download was another
# chance for a flaky upstream to fail. Generated files -- trimmed chain segments, rendered
# HTML, TM-align output -- stay in the per-session workdir, because they are that
# session's output rather than a copy of something upstream.

_STRUCTURE_DIR: Path | None = None
_STRUCTURE_DIR_LOCK = threading.Lock()
_FILE_LOCKS: dict[str, threading.Lock] = {}


def shared_structure_dir() -> Path:
    """A process-wide directory for downloaded structure files."""
    global _STRUCTURE_DIR
    with _STRUCTURE_DIR_LOCK:
        if _STRUCTURE_DIR is None:
            base = os.getenv("SCOP3P_STRUCTURE_CACHE_DIR")
            _STRUCTURE_DIR = (
                Path(base) if base else Path(tempfile.gettempdir()) / "scop3p_structures"
            )
            _STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
            LOGGER.info(
                "structure cache directory=%s",
                _STRUCTURE_DIR,
                extra={"event": "cache"},
            )
        return _STRUCTURE_DIR


def structure_file_lock(filename: str) -> threading.Lock:
    """One lock per filename, so two sessions do not download the same file at once.

    Without it, two tabs asking for 2IVT simultaneously both fetch it and both write the
    same path -- wasted bandwidth, and a reader can see a half-written file.
    """
    with _STRUCTURE_DIR_LOCK:
        return _FILE_LOCKS.setdefault(filename, threading.Lock())


def cached_structure_file(filename: str, minimum_bytes: int = 1) -> Path | None:
    """An already-downloaded file, if it is present and non-trivial."""
    path = shared_structure_dir() / filename
    if path.exists() and path.stat().st_size >= minimum_bytes:
        return path
    return None
