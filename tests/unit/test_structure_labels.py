"""One vocabulary for describing a selectable structure, across every protocol.

Three apps let a user pick a structure and each had grown its own notation for the same
facts. The point of this module is that a user moving between protocols never has to
re-read the notation, so the tests that matter most are the cross-app consistency ones at
the bottom: they fail if any app goes back to formatting its own labels.
"""

from __future__ import annotations

import pytest

from common.structure_labels import (
    ALL_CHAINS_PLACEHOLDER,
    ALPHAFOLD_OPTION_LABEL,
    CHOOSE_ENTRY_PLACEHOLDER,
    NO_PROTEIN_PLACEHOLDER,
    NO_STRUCTURES_PLACEHOLDER,
    SEPARATOR,
    chain_label,
    chain_option_label,
    format_chains,
    format_range,
    format_resolution,
    structure_option_label,
)


# --------------------------------------------------------------------------------------
# Resolution, which arrives in every shape the upstreams can manage
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value, expected",
    [
        (2.0, "2.00 A"),          # PDBe returns a float
        (2, "2.00 A"),
        ("2.00 A", "2.00 A"),     # UniProt returns a formatted string
        ("2.0A", "2.00 A"),
        ("2", "2.00 A"),
        (None, ""),
        ("", ""),
        ("-", ""),                # methods with no resolution report a dash
        ("?", ""),
        ("nan", ""),
    ],
)
def test_resolution_is_normalised(value, expected) -> None:
    assert format_resolution(value) == expected


def test_an_unparseable_resolution_is_passed_through_not_dropped() -> None:
    """Better to show something odd than to silently hide a value the upstream sent."""
    assert format_resolution("low-resolution EM") == "low-resolution EM"


# --------------------------------------------------------------------------------------
# Chains
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value, expected",
    [
        ({"A": None}, "chain A"),
        ({"A": None, "B": None}, "chains A, B"),
        ({"B": None, "A": None}, "chains A, B"),   # sorted, so labels are stable
        (["A", "B"], "chains A, B"),
        ("A", "chain A"),
        ({}, ""),
        ([], ""),
        (None, ""),
    ],
)
def test_chains_are_described_and_sorted(value, expected) -> None:
    assert format_chains(value) == expected


def test_singular_and_plural_are_not_the_same_word() -> None:
    """"chains A" reads like a defect."""
    assert format_chains({"A": None}).startswith("chain ")
    assert format_chains({"A": None, "B": None}).startswith("chains ")


# --------------------------------------------------------------------------------------
# Ranges
# --------------------------------------------------------------------------------------


def test_a_range_is_labelled_with_its_coordinate_system() -> None:
    assert format_range(705, 1013) == "UniProt 705-1013"


@pytest.mark.parametrize("start, end", [("?", "?"), ("?", 100), (1, "?"), (None, None)])
def test_an_unknown_bound_produces_nothing(start, end) -> None:
    """The upstreams use "?" for an unknown bound.

    Rendering it would put "UniProt ?-1013" in front of the user, which says nothing and
    looks broken.
    """
    assert format_range(start, end) == ""


# --------------------------------------------------------------------------------------
# Entry-level labels
# --------------------------------------------------------------------------------------


def test_an_entry_label_reads_id_method_resolution_chains() -> None:
    label = structure_option_label(
        "2ivs", method="X-ray", resolution="2.00 A", chains={"A": None, "B": None}
    )
    assert label == "2IVS · X-ray · 2.00 A · chains A, B"


def test_the_pdb_id_is_always_first_and_upper_case() -> None:
    assert structure_option_label("2ivs").startswith("2IVS")


def test_an_entry_with_nothing_but_an_id_is_still_usable() -> None:
    assert structure_option_label("9zzz") == "9ZZZ"


def test_missing_parts_leave_no_dangling_separator() -> None:
    for label in (
        structure_option_label("9ZZZ"),
        structure_option_label("9ZZZ", method="", resolution=None, chains=None),
        structure_option_label("9ZZZ", resolution="-"),
        chain_option_label("9ZZZ", "A", unp_start="?", unp_end="?"),
    ):
        assert not label.endswith(SEPARATOR.rstrip())
        assert not label.startswith(SEPARATOR.strip())
        assert f"{SEPARATOR}{SEPARATOR}" not in label
        assert "··" not in label


def test_coverage_is_rendered_as_a_percentage() -> None:
    label = structure_option_label("2IVS", coverage=0.29)
    assert "29% cover" in label


@pytest.mark.parametrize("coverage", ["", "n/a", None])
def test_an_unusable_coverage_is_omitted(coverage) -> None:
    assert "cover" not in structure_option_label("2IVS", coverage=coverage)


# --------------------------------------------------------------------------------------
# Chain-level labels
# --------------------------------------------------------------------------------------


def test_a_chain_level_label_names_the_entry_and_the_chain() -> None:
    label = chain_option_label(
        "2ivs", "A", unp_start=705, unp_end=1013, method="X-ray", resolution=2.0
    )
    assert label == "2IVS · chain A · UniProt 705-1013 · X-ray · 2.00 A"


def test_a_chain_picker_label_carries_the_range() -> None:
    """The range is what shows a structure covers a domain, not the whole protein."""
    assert chain_label("A", (705, 1013)) == "A (705-1013)"


def test_a_chain_with_no_known_range_is_just_the_chain() -> None:
    assert chain_label("A", None) == "A"
    assert chain_label("A", ()) == "A"


def test_a_malformed_range_does_not_raise() -> None:
    assert chain_label("A", ("x", "y")) == "A"


# --------------------------------------------------------------------------------------
# Cross-app consistency -- the reason this module exists
# --------------------------------------------------------------------------------------


def test_every_protocol_describes_the_same_entry_identically() -> None:
    """structure_viz and topology_viewer must not diverge again.

    They reach the formatter from different upstreams -- UniProt cross-references give a
    resolution string, PDBe gives a float -- so this also pins that both normalise.
    """
    from_uniprot = structure_option_label(
        "2IVS", method="X-ray", resolution="2.00 A", chains={"A": (705, 1013), "B": (705, 1013)}
    )
    from_pdbe = structure_option_label(
        "2ivs", method="X-ray", resolution=2.0, chains={"B": (705, 1013), "A": (705, 1013)}
    )
    assert from_uniprot == from_pdbe == "2IVS · X-ray · 2.00 A · chains A, B"


def test_the_entry_and_chain_labels_share_one_notation() -> None:
    """The two granularities must read as one family, not two conventions."""
    entry = structure_option_label("2IVS", method="X-ray", resolution=2.0, chains={"A": None})
    chain = chain_option_label("2IVS", "A", method="X-ray", resolution=2.0)
    for label in (entry, chain):
        assert label.startswith("2IVS")
        assert SEPARATOR in label
        assert "2.00 A" in label
        assert "chain A" in label


def test_the_alphafold_option_has_one_name() -> None:
    """Each app naming this differently was the most visible inconsistency."""
    assert ALPHAFOLD_OPTION_LABEL
    assert "AlphaFold" in ALPHAFOLD_OPTION_LABEL


def test_placeholders_are_distinct_and_explain_themselves() -> None:
    placeholders = [
        NO_PROTEIN_PLACEHOLDER,
        NO_STRUCTURES_PLACEHOLDER,
        CHOOSE_ENTRY_PLACEHOLDER,
        ALL_CHAINS_PLACEHOLDER,
    ]
    assert len(set(placeholders)) == len(placeholders)
    for text in placeholders:
        assert text.strip(), "an empty placeholder renders as a blank option"


def test_no_structures_is_not_the_same_message_as_no_protein() -> None:
    """"Set a protein first" and "this protein has no PDB entries" are different facts."""
    assert NO_PROTEIN_PLACEHOLDER != NO_STRUCTURES_PLACEHOLDER
