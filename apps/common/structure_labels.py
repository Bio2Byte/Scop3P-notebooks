"""One way to describe a selectable structure, shared by every protocol.

Three apps let the user pick a structure, and each had grown its own wording for the
same facts:

    structure_viz     2IVS · X-ray · 2.00 A · chains A, B
    topology_viewer   2IVS · X-ray · 2.00 A · 29% cover
    rinalign          2IVS chain A [705-1013] (2.0A)

Same protein, same entry, three formats -- and rinalign built its own in three separate
code paths. A user moving between protocols had to re-learn the notation each time, and
"is this the same 2IVS I was looking at?" is a question the interface should never raise.

This module is the single source of that wording. It is deliberately dependency-free --
no shiny, no pandas -- so the service layer and the UI layer can both use it.

Two levels of granularity, because the protocols genuinely differ:

* Entry level (`structure_option_label`): the user picks an entry, then a chain.
  structure_viz and topology_viewer work this way.
* Chain level (`chain_option_label`): the user picks an entry-and-chain in one go,
  because a residue interaction network is built on exactly one chain. rinalign works
  this way.

Both start with the PDB id and separate facts with the same middle dot, so the two read
as members of one family rather than as unrelated notations.
"""

from __future__ import annotations

SEPARATOR = " · "

#: What every protocol calls the AlphaFold option. It is not a PDB entry, so it has no
#: method or resolution, and each app naming it differently was the most visible part of
#: the inconsistency.
ALPHAFOLD_OPTION_LABEL = "AlphaFold model (full length)"

#: Shown by an entry picker before an accession has been set.
NO_PROTEIN_PLACEHOLDER = "Set a protein first"

#: Shown by an entry picker when the accession has no cross-referenced structures.
NO_STRUCTURES_PLACEHOLDER = "No PDB entries for this protein"

#: The empty option that invites a choice, rather than silently pre-selecting one.
CHOOSE_ENTRY_PLACEHOLDER = "-- select a PDB entry --"

#: The "no particular chain" option, where a protocol can work on the whole entry.
ALL_CHAINS_PLACEHOLDER = "All chains"

#: Shown by a chain picker whose options only exist once structures have been loaded. An
#: empty select renders as a blank box, which reads as a broken control rather than as
#: "nothing to choose yet".
NO_STRUCTURES_LOADED_PLACEHOLDER = "Load structures first"

#: Shown when the lookup that fills a picker *failed*, as opposed to succeeding with
#: nothing. These must never read the same: "no entries" is a fact about the protein and
#: invites the user to move on, while a failed lookup is a broken run that a retry may
#: well fix. The upstreams here fail intermittently, so this is a common case, and the
#: text has to say what to do about it because the picker is the only place the user
#: looks.
LOOKUP_FAILED_PLACEHOLDER = "Lookup failed - press Set protein to retry"

#: The same, for a picker filled by a button labelled "Fetch".
LOOKUP_FAILED_RETRY_FETCH = "Lookup failed - press Fetch to retry"


def format_resolution(resolution: object) -> str:
    """Resolution as ``2.00 A``.

    Accepts what the various upstreams actually return: a float from PDBe, a string like
    ``"2.00 A"`` from UniProt's cross-reference properties, ``"-"`` for methods that have
    no resolution, or nothing at all.
    """
    if resolution is None:
        return ""
    if isinstance(resolution, (int, float)):
        return f"{float(resolution):.2f} A"

    text = str(resolution).strip()
    if not text or text in {"-", "?", "None", "nan"}:
        return ""
    # "2.00 A", "2.0A" and "2" all arrive; normalise the ones that parse.
    numeric = text.rstrip("Aa ").strip()
    try:
        return f"{float(numeric):.2f} A"
    except ValueError:
        return text


def _joined(pdb_id: str, *parts: str) -> str:
    bits = [str(pdb_id).upper().strip()]
    bits.extend(part for part in parts if part)
    return SEPARATOR.join(bits)


def format_chains(chains: object) -> str:
    """``chains A, B`` from a mapping, a sequence, or a single chain id."""
    if not chains:
        return ""
    if isinstance(chains, str):
        names = [chains]
    elif isinstance(chains, dict):
        names = sorted(str(name) for name in chains)
    else:
        names = sorted(str(name) for name in chains)
    names = [name for name in names if name]
    if not names:
        return ""
    return f"chains {', '.join(names)}" if len(names) > 1 else f"chain {names[0]}"


def format_range(start: object, end: object) -> str:
    """``UniProt 705-1013``, or nothing when either bound is unknown.

    The upstreams use ``"?"`` for an unknown bound, which must not be rendered: a label
    reading "UniProt ?-1013" tells the user nothing and looks like a defect.
    """
    try:
        return f"UniProt {int(start)}-{int(end)}"
    except (TypeError, ValueError):
        return ""


def structure_option_label(
    pdb_id: str,
    *,
    method: object = "",
    resolution: object = None,
    chains: object = None,
    coverage: object = None,
) -> str:
    """Describe one PDB entry for an entry-level picker.

    ``2IVS · X-ray · 2.00 A · chains A, B``

    Every part is optional, because the upstreams are inconsistent about what they
    report; an entry with nothing but an id still gets a usable label.
    """
    coverage_text = ""
    if coverage is not None:
        try:
            coverage_text = f"{float(coverage) * 100:.0f}% cover"
        except (TypeError, ValueError):
            coverage_text = ""
    return _joined(
        pdb_id,
        str(method or "").strip(),
        format_resolution(resolution),
        format_chains(chains),
        coverage_text,
    )


def chain_option_label(
    pdb_id: str,
    chain: str,
    *,
    unp_start: object = None,
    unp_end: object = None,
    method: object = "",
    resolution: object = None,
) -> str:
    """Describe one chain of one entry, for a picker that selects both at once.

    ``2IVS · chain A · UniProt 705-1013 · X-ray · 2.00 A``
    """
    return _joined(
        pdb_id,
        format_chains(chain),
        format_range(unp_start, unp_end),
        str(method or "").strip(),
        format_resolution(resolution),
    )


def chain_label(chain: str, uniprot_range: tuple[int, int] | None = None) -> str:
    """A chain on its own, for a chain picker that follows an entry picker.

    ``A (705-1013)``. The range is what tells the user that a structure covers one
    domain rather than the whole protein.
    """
    name = str(chain).strip()
    if not uniprot_range:
        return name
    try:
        return f"{name} ({int(uniprot_range[0])}-{int(uniprot_range[1])})"
    except (TypeError, ValueError, IndexError):
        return name
