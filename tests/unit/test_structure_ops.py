from __future__ import annotations

from pathlib import Path
import pytest

from common.structure_viz import StructureVizService
from common.structure_viz import StructureOps


PDB_MINI = """ATOM      1  N   ALA A   1      11.104  13.207   2.100  1.00 10.00           N
ATOM      2  CA  ALA A   1      12.200  12.300   2.300  1.00 10.00           C
ATOM      3  N   GLY A   2      13.104  14.207   3.100  1.00 10.00           N
ATOM      4  CA  GLY A   2      14.200  14.300   3.300  1.00 10.00           C
TER
END
"""

PDB_MULTI_CHAIN = """ATOM      1  N   ALA A   5      11.104  13.207   2.100  1.00 10.00           N
ATOM      2  CA  ALA A   5      12.200  12.300   2.300  1.00 10.00           C
ATOM      3  N   GLY A  10      13.104  14.207   3.100  1.00 10.00           N
ATOM      4  CA  GLY A  10      14.200  14.300   3.300  1.00 10.00           C
TER
ATOM      5  N   SER B 200      21.104  23.207   2.100  1.00 10.00           N
ATOM      6  CA  SER B 200      22.200  22.300   2.300  1.00 10.00           C
ATOM      7  N   TYR B 210      23.104  24.207   3.100  1.00 10.00           N
ATOM      8  CA  TYR B 210      24.200  24.300   3.300  1.00 10.00           C
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


def test_chain_ranges_from_pdb_returns_all_available_ranges(tmp_path: Path) -> None:
    source = tmp_path / "multi_chain.pdb"
    source.write_text(PDB_MULTI_CHAIN)

    ranges = StructureOps.chain_ranges_from_pdb(source)
    assert ranges == {"A": (5, 10), "B": (200, 210)}


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


def test_run_tmalign_raises_runtime_error_with_subprocess_output(monkeypatch, tmp_path: Path) -> None:
    pdb1 = tmp_path / "a.pdb"
    pdb2 = tmp_path / "b.pdb"
    pdb1.write_text(PDB_MINI)
    pdb2.write_text(PDB_MINI)

    def _mock_run(cmd, cwd, capture_output, text, check):  # noqa: ANN001
        error = __import__("subprocess").CalledProcessError(139, cmd, output="")
        error.stderr = "segmentation fault"
        raise error

    monkeypatch.setattr("apps.common.structure_viz.subprocess.run", _mock_run)

    with pytest.raises(RuntimeError, match="segmentation fault"):
        StructureOps.run_tmalign(pdb1, pdb2, tmp_path, out_name="aligned")


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


def test_validate_pdb_id_accepts_2ivt_and_rejects_bad_id() -> None:
    assert StructureOps.validate_pdb_id("2IVT") == "2IVT"
    with pytest.raises(ValueError, match="Expected a 4-character RCSB identifier"):
        StructureOps.validate_pdb_id("21VTX")


def test_save_chain_segment_rejects_missing_chain(tmp_path: Path) -> None:
    source = tmp_path / "source.pdb"
    source.write_text(PDB_MINI)

    with pytest.raises(ValueError, match="Chain 'B' not found"):
        StructureOps.save_chain_segment(source, tmp_path / "seg.pdb", "B", 1, 2)


def test_save_chain_segment_rejects_empty_range(tmp_path: Path) -> None:
    source = tmp_path / "source.pdb"
    source.write_text(PDB_MINI)

    with pytest.raises(ValueError, match="outside chain 'A' range 1-2"):
        StructureOps.save_chain_segment(source, tmp_path / "seg.pdb", "A", 5, 10)


def test_download_pdb_uses_uppercase_2ivt_url(monkeypatch, tmp_path: Path) -> None:
    service = StructureVizService(tmp_path)

    class _Response:
        content = PDB_MINI.encode()

        def raise_for_status(self) -> None:
            return None

    def _mock_get(url, timeout):  # noqa: ANN001
        assert url == "https://files.rcsb.org/download/2IVT.pdb"
        return _Response()

    monkeypatch.setattr("apps.common.structure_viz.requests.get", _mock_get)

    downloaded = service.download_pdb("2ivt")
    assert downloaded.name == "2IVT.pdb"
    assert downloaded.read_text() == PDB_MINI
