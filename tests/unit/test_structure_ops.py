from __future__ import annotations

from pathlib import Path

from common.structure_viz import StructureVizService
from common.structure_viz import StructureOps


PDB_MINI = """ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00 10.00           N
ATOM      2  CA  ALA A   1      12.200  12.300   2.300  1.00 10.00           C
ATOM      3  N   GLY A   2      13.104  14.207   3.100  1.00 10.00           N
ATOM      4  CA  GLY A   2      14.200  14.300   3.300  1.00 10.00           C
TER
END
"""


def test_bfactor_pdb_rewrites_selected_metric(tmp_path: Path) -> None:
    source = tmp_path / "in.pdb"
    source.write_text(PDB_MINI)
    out = tmp_path / "out.pdb"

    import pandas as pd

    df = pd.DataFrame({"metric": [1.23, 9.87]})
    result = StructureOps.bfactor_pdb(source, df, "metric", out, chain="A")
    lines = result.read_text().splitlines()
    assert any("  1.23" in line for line in lines if line.startswith("ATOM"))
    assert any("  9.87" in line for line in lines if line.startswith("ATOM"))


def test_chain_range_from_pdb(tmp_path: Path) -> None:
    source = tmp_path / "range.pdb"
    source.write_text(PDB_MINI)
    start, end = StructureOps.chain_range_from_pdb(source, "A")
    assert start == 1
    assert end == 2


def test_run_tmalign_finds_existing_candidate(monkeypatch, tmp_path: Path) -> None:
    pdb1 = tmp_path / "a.pdb"
    pdb2 = tmp_path / "b.pdb"
    pdb1.write_text(PDB_MINI)
    pdb2.write_text(PDB_MINI)
    aligned = tmp_path / "aligned"
    aligned.write_text("MODEL\nEND\n")

    class _Result:
        stdout = "Aligned\n"

    def _mock_run(cmd, cwd, capture_output, text, check):  # noqa: ANN001
        assert cwd == tmp_path
        return _Result()

    monkeypatch.setattr("apps.common.structure_viz.subprocess.run", _mock_run)

    path, report = StructureOps.run_tmalign(pdb1, pdb2, tmp_path, out_name="aligned")
    assert path == aligned
    assert report == "Aligned\n"


def test_resolve_uploaded_or_remote_pdb_prefers_upload(tmp_path: Path) -> None:
    service = StructureVizService(tmp_path)
    source = tmp_path / "upload_source.pdb"
    source.write_text(PDB_MINI)

    resolved = service.resolve_uploaded_or_remote_pdb(
        [{"datapath": str(source), "name": "user_file.pdb"}],
        "2IVT",
        target_name="copied.pdb",
    )

    assert resolved == tmp_path / "copied.pdb"
    assert resolved.read_text() == PDB_MINI
