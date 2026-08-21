"""The process-wide lookup cache.

The properties pinned here are the ones whose absence would be invisible: a cache that
looks like it works while sharing nothing, or one that turns a blip into a dead protocol.
Two of these tests exist because the bug happened during development, not in theory.
"""

from __future__ import annotations

import threading
import time

import pandas as pd
import pytest

from common.cache import (
    REGISTRY,
    LookupCache,
    cache_report,
    clear_all,
    memoize,
    shared_structure_dir,
    structure_file_lock,
)


# --------------------------------------------------------------------------------------
# The store
# --------------------------------------------------------------------------------------


def test_a_value_round_trips() -> None:
    cache = LookupCache(name="t")
    cache.set("k", 42)
    assert cache.get("k") == (True, 42)


def test_a_miss_is_distinguishable_from_a_cached_none() -> None:
    """The reason get() returns a tuple rather than a sentinel."""
    cache = LookupCache(name="t")
    assert cache.get("absent") == (False, None)
    cache.set("present", None)
    assert cache.get("present") == (True, None)


def test_an_entry_expires() -> None:
    cache = LookupCache(name="t", ttl_seconds=0.05)
    cache.set("k", 1)
    time.sleep(0.08)
    found, _ = cache.get("k")
    assert not found
    assert cache.stats.expirations == 1


def test_the_least_recently_used_entry_is_evicted() -> None:
    cache = LookupCache(name="t", max_entries=2)
    cache.set("a", 1)
    cache.set("b", 2)
    cache.get("a")           # "a" is now the most recent
    cache.set("c", 3)        # evicts "b"
    assert cache.get("a")[0]
    assert cache.get("c")[0]
    assert not cache.get("b")[0]
    assert cache.stats.evictions == 1


def test_the_bound_is_actually_enforced() -> None:
    cache = LookupCache(name="t", max_entries=10)
    for index in range(100):
        cache.set(index, index)
    assert len(cache) == 10


def test_invalidate_can_target_a_subset() -> None:
    cache = LookupCache(name="t")
    cache.set(("P1",), 1)
    cache.set(("P2",), 2)
    removed = cache.invalidate(lambda key: key == ("P1",))
    assert removed == 1
    assert not cache.get(("P1",))[0]
    assert cache.get(("P2",))[0]


def test_concurrent_access_does_not_corrupt_the_store() -> None:
    """Handlers run under asyncio.to_thread, so this is a real access pattern."""
    cache = LookupCache(name="t", max_entries=1000)

    def worker(offset: int) -> None:
        for index in range(200):
            cache.set((offset, index), index)
            cache.get((offset, index))

    threads = [threading.Thread(target=worker, args=(n,)) for n in range(8)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert len(cache) <= 1000


# --------------------------------------------------------------------------------------
# The decorator
# --------------------------------------------------------------------------------------


def test_a_repeated_call_does_not_reach_the_function() -> None:
    calls = {"n": 0}

    @memoize(name="test.basic")
    def fetch(accession: str) -> str:
        calls["n"] += 1
        return f"seq-{accession}"

    assert fetch("P1") == "seq-P1"
    assert fetch("P1") == "seq-P1"
    assert calls["n"] == 1


def test_a_failure_is_never_cached() -> None:
    """The property that matters most.

    These upstreams fail transiently. Caching the exception would make the retry button
    useless and the only cure a restart, turning a blip into an outage.
    """
    calls = {"n": 0}

    @memoize(name="test.failure")
    def flaky(accession: str) -> str:
        calls["n"] += 1
        if calls["n"] == 1:
            raise RuntimeError("upstream dropped the connection")
        return "recovered"

    with pytest.raises(RuntimeError):
        flaky("P1")
    # The retry must actually reach the function.
    assert flaky("P1") == "recovered"
    assert calls["n"] == 2
    assert REGISTRY["test.failure"].stats.errors_not_cached == 1


def test_the_key_ignores_self_so_sessions_share() -> None:
    """Each browser session builds its own service instance.

    Keying on the instance would give every session a private cache: no sharing at all,
    while still looking like a working cache.
    """
    calls = {"n": 0}

    class Service:
        @memoize(name="test.skipself")
        def fetch(self, accession: str) -> str:
            calls["n"] += 1
            return "value"

    Service().fetch("P1")
    Service().fetch("P1")  # a different instance, as a second session would be
    assert calls["n"] == 1


def test_two_functions_sharing_a_name_share_one_store() -> None:
    """Mutation Effect and Structure Visualisation fetch the same UniProt FASTA.

    This regressed during development: the decorator built a fresh cache per function, so
    the shared name shared nothing and orphaned all but the last from the registry -- which
    also meant clear_all() silently missed it.
    """
    calls = {"a": 0, "b": 0}

    @memoize(name="test.shared")
    def first(accession: str) -> str:
        calls["a"] += 1
        return "seq"

    @memoize(name="test.shared")
    def second(accession: str) -> str:
        calls["b"] += 1
        return "seq"

    first("P1")
    second("P1")
    assert calls == {"a": 1, "b": 0}, "the second function did not see the first's entry"
    assert first.cache is second.cache


def test_every_cache_is_reachable_from_the_registry() -> None:
    """clear_all() works through REGISTRY, so an unregistered cache is never cleared."""

    @memoize(name="test.registered")
    def fetch(accession: str) -> str:
        return "value"

    fetch("P1")
    assert REGISTRY["test.registered"] is fetch.cache
    assert len(fetch.cache) == 1
    clear_all()
    assert len(fetch.cache) == 0


def test_arguments_are_part_of_the_key() -> None:
    @memoize(name="test.keys")
    def fetch(accession: str) -> str:
        return f"seq-{accession}"

    assert fetch("P1") == "seq-P1"
    assert fetch("P2") == "seq-P2"


# --------------------------------------------------------------------------------------
# Copy-on-return
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("call_index", [1, 2, 3])
def test_a_cached_dataframe_cannot_be_mutated_by_any_caller(call_index: int) -> None:
    """Two sessions hold the same object otherwise, and one can corrupt the other.

    Parametrised over which call does the mutating on purpose. Copying only on the way
    into the cache protects the *first* caller's result and leaves every later caller
    holding the live object -- so a test that only mutates the first result passes while
    the hit path is wide open. That is exactly what happened here.
    """

    @memoize(name="test.frame")
    def fetch(accession: str) -> pd.DataFrame:
        return pd.DataFrame({"position": [1, 2]})

    for _ in range(call_index - 1):
        fetch("P1")
    fetch("P1").loc[0, "position"] = 999
    assert fetch("P1").loc[0, "position"] == 1


@pytest.mark.parametrize("call_index", [1, 2, 3])
def test_a_cached_list_cannot_be_reordered_by_any_caller(call_index: int) -> None:
    """A caller sorting the returned list would reorder it for everyone else."""

    @memoize(name="test.list")
    def fetch(accession: str) -> list[str]:
        return ["b", "a", "c"]

    for _ in range(call_index - 1):
        fetch("P1")
    fetch("P1").sort()
    assert fetch("P1") == ["b", "a", "c"]


def test_the_object_handed_out_is_never_the_stored_one() -> None:
    """Directly, rather than through a mutation, so the guarantee is unambiguous."""

    @memoize(name="test.identity")
    def fetch(accession: str) -> list[str]:
        return ["a"]

    first, second = fetch("P1"), fetch("P1")
    assert first is not second
    assert first == second


@pytest.mark.parametrize(
    "value", [{"a": 1}, {"a", "b"}, pd.Series([1, 2])]
)
def test_other_mutable_results_are_copied(value) -> None:
    @memoize(name="test.mutable")
    def fetch(accession: str):
        return value

    assert fetch("P1") is not value


def test_an_immutable_result_is_returned_as_is() -> None:
    """Copying a string per call would be pure waste."""

    @memoize(name="test.str")
    def fetch(accession: str) -> str:
        return "MVLSPADKTN"

    assert fetch("P1") == "MVLSPADKTN"


# --------------------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------------------


def test_hits_and_misses_are_counted() -> None:
    @memoize(name="test.stats")
    def fetch(accession: str) -> str:
        return "value"

    fetch("P1")
    fetch("P1")
    report = cache_report()["test.stats"]
    assert report["hits"] == 1
    assert report["misses"] == 1


def test_the_report_covers_the_real_caches() -> None:
    """A sanity check that the protocols' lookups are registered at all."""
    import common.mutation_effect  # noqa: F401
    import common.rinalign  # noqa: F401
    import common.services  # noqa: F401
    import common.structure_viz  # noqa: F401

    names = set(cache_report())
    for expected in (
        "uniprot.sequence.fasta",
        "uniprot.pdb.xrefs",
        "uniprot.disease.variants",
        "uniprot.ptm.features",
        "scop3p.modifications",
    ):
        assert expected in names, f"{expected} is not memoized"


# --------------------------------------------------------------------------------------
# Structure files
# --------------------------------------------------------------------------------------


def test_the_structure_directory_is_shared_and_writable() -> None:
    directory = shared_structure_dir()
    assert directory.exists()
    assert shared_structure_dir() is directory or shared_structure_dir() == directory


def test_the_structure_directory_honours_its_environment_override(monkeypatch, tmp_path) -> None:
    import common.cache as cache_module

    monkeypatch.setenv("SCOP3P_STRUCTURE_CACHE_DIR", str(tmp_path / "elsewhere"))
    monkeypatch.setattr(cache_module, "_STRUCTURE_DIR", None, raising=False)
    assert shared_structure_dir() == tmp_path / "elsewhere"


def test_one_lock_per_filename() -> None:
    """Two sessions asking for the same entry must serialise on it, not on everything."""
    assert structure_file_lock("2IVT.pdb") is structure_file_lock("2IVT.pdb")
    assert structure_file_lock("2IVT.pdb") is not structure_file_lock("1A3N.pdb")


# --------------------------------------------------------------------------------------
# Bio2Byte predictions
# --------------------------------------------------------------------------------------
# The most expensive thing the toolkit computes -- about 16 seconds for 1100 residues --
# and deterministic, so worth caching. It also carries the sharpest correctness risk in
# the whole cache, because Mutation Effect predicts a *mutant* under the wild type's
# accession.


WILD_TYPE = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQ"
MUTANT = WILD_TYPE[:10] + "W" + WILD_TYPE[11:]


def _service(monkeypatch):
    """A MutationEffectService whose predictor is a counter, not b2bTools."""
    from common.mutation_effect import MutationEffectService

    calls: list[str] = []

    def fake_predict(self, accession: str, sequence: str) -> dict:
        calls.append(sequence)
        return {"proteins": {accession: {"seq": sequence, "n": len(calls)}}}

    monkeypatch.setattr(
        MutationEffectService, "_predict_biophysical_uncached", fake_predict
    )
    return MutationEffectService(), calls


def test_the_wild_type_prediction_is_cached(monkeypatch) -> None:
    service, calls = _service(monkeypatch)
    service.predict_biophysical("P1", WILD_TYPE)
    service.predict_biophysical("P1", WILD_TYPE)
    assert calls == [WILD_TYPE], "the wild type was predicted twice"


def test_a_mutant_prediction_is_not_cached(monkeypatch) -> None:
    """A mutant is seen once, so a slot spent on it never earns a hit."""
    service, calls = _service(monkeypatch)
    service.predict_biophysical("P1", MUTANT, wild_type=False)
    service.predict_biophysical("P1", MUTANT, wild_type=False)
    assert calls == [MUTANT, MUTANT], "a mutant prediction was cached"


def test_a_mutant_never_receives_the_wild_type_prediction(monkeypatch) -> None:
    """The correctness risk that makes this cache worth testing carefully.

    Same accession, different sequence. Keying on the accession alone would return the
    wild type's numbers for the mutant -- silently destroying the very comparison the
    Mutation Effect protocol exists to make, with no error anywhere.
    """
    service, calls = _service(monkeypatch)
    wild = service.predict_biophysical("P1", WILD_TYPE)
    mutant = service.predict_biophysical("P1", MUTANT, wild_type=False)

    assert wild["proteins"]["P1"]["seq"] == WILD_TYPE
    assert mutant["proteins"]["P1"]["seq"] == MUTANT
    assert calls == [WILD_TYPE, MUTANT]


def test_predicting_mutants_does_not_evict_the_wild_type(monkeypatch) -> None:
    """Exploring mutations must not cost the one entry every comparison needs."""
    service, calls = _service(monkeypatch)
    service.predict_biophysical("P1", WILD_TYPE)
    for index in range(50):  # more than the prediction cache bound
        variant = WILD_TYPE[:index + 1] + "W" + WILD_TYPE[index + 2:]
        service.predict_biophysical("P1", variant, wild_type=False)

    calls.clear()
    service.predict_biophysical("P1", WILD_TYPE)
    assert calls == [], "the wild-type prediction was evicted by mutant traffic"


def test_a_changed_canonical_sequence_is_not_served_stale(monkeypatch) -> None:
    """The sequence is in the key, so an upstream revision re-predicts."""
    service, calls = _service(monkeypatch)
    service.predict_biophysical("P1", WILD_TYPE)
    service.predict_biophysical("P1", WILD_TYPE + "AAA")
    assert len(calls) == 2


def test_a_failed_prediction_is_not_cached(monkeypatch) -> None:
    """b2bTools can fail; the button must remain a retry."""
    from common.mutation_effect import MutationEffectService

    attempts = {"n": 0}

    def flaky(self, accession: str, sequence: str) -> dict:
        attempts["n"] += 1
        if attempts["n"] == 1:
            raise RuntimeError("predictor exploded")
        return {"proteins": {accession: {"seq": sequence}}}

    monkeypatch.setattr(MutationEffectService, "_predict_biophysical_uncached", flaky)
    service = MutationEffectService()
    with pytest.raises(RuntimeError):
        service.predict_biophysical("P1", WILD_TYPE)
    assert service.predict_biophysical("P1", WILD_TYPE)["proteins"]["P1"]["seq"] == WILD_TYPE
    assert attempts["n"] == 2


def test_the_two_prediction_caches_stay_separate() -> None:
    """One returns a DataFrame, the other a raw dict.

    Sharing a name would hand a caller the other shape -- an AttributeError deep in a
    render function, a long way from here.
    """
    import common.mutation_effect  # noqa: F401
    import common.structure_viz  # noqa: F401

    names = set(cache_report())
    assert {"b2b.prediction.table", "b2b.prediction.raw"} <= names
    assert REGISTRY["b2b.prediction.table"] is not REGISTRY["b2b.prediction.raw"]


def test_prediction_caches_are_bounded_tighter_than_json_lookups() -> None:
    """A prediction frame is orders of magnitude larger than a JSON payload."""
    from common.cache import DEFAULT_MAX_ENTRIES, PREDICTION_MAX_ENTRIES

    assert PREDICTION_MAX_ENTRIES < DEFAULT_MAX_ENTRIES
    assert REGISTRY["b2b.prediction.raw"].max_entries == PREDICTION_MAX_ENTRIES


def test_a_mutant_is_safe_even_if_the_caller_forgets_to_opt_out(monkeypatch) -> None:
    """Defence in depth for the same correctness risk.

    The bypass keeps mutants out of the cache, but the bypass is a call-site decision and
    call sites get edited. If someone predicts a mutant on the cached path by mistake, the
    sequence being part of the key is what still prevents the wild type's numbers being
    returned for it. Both protections are wanted, because only this one survives a caller
    getting it wrong.
    """
    service, calls = _service(monkeypatch)
    wild = service.predict_biophysical("P1", WILD_TYPE)
    # Deliberately on the cached path, as a mistaken call site would be.
    mutant = service.predict_biophysical("P1", MUTANT)

    assert wild["proteins"]["P1"]["seq"] == WILD_TYPE
    assert mutant["proteins"]["P1"]["seq"] == MUTANT
    assert calls == [WILD_TYPE, MUTANT], "the mutant was served the wild type's prediction"
