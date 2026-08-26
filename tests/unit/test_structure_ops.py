from __future__ import annotations

from pathlib import Path
import pytest

from common.structure_viz import StructureVizService
from common.structure_viz import StructureOps
from structure_viz.app import _tm_source_signature


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

    monkeypatch.setattr("common.structure_viz.subprocess.run", _mock_run)

    output = StructureOps.run_tmalign(pdb1, pdb2, tmp_path, out_name="aligned")
    assert output.superposed == aligned
    assert output.report == "Aligned\n"
    # The reference half of the superposition is the second input, unmoved. TM-align never
    # writes the two structures as one file, so both paths have to come back or the viewer
    # can only ever draw one of them.
    assert output.reference == pdb2


def test_run_tmalign_raises_runtime_error_with_subprocess_output(monkeypatch, tmp_path: Path) -> None:
    pdb1 = tmp_path / "a.pdb"
    pdb2 = tmp_path / "b.pdb"
    pdb1.write_text(PDB_MINI)
    pdb2.write_text(PDB_MINI)

    def _mock_run(cmd, cwd, capture_output, text, check):  # noqa: ANN001
        error = __import__("subprocess").CalledProcessError(139, cmd, output="")
        error.stderr = "segmentation fault"
        raise error

    monkeypatch.setattr("common.structure_viz.subprocess.run", _mock_run)

    with pytest.raises(RuntimeError, match="segmentation fault"):
        StructureOps.run_tmalign(pdb1, pdb2, tmp_path, out_name="aligned")


MATRIX_IDENTITY_SHIFT = """ ------ The rotation matrix to rotate Chain_1 to Chain_2 ------
m               t[m]        u[m][0]        u[m][1]        u[m][2]
0     1.0000000000   1.0000000000   0.0000000000   0.0000000000
1     2.0000000000   0.0000000000   1.0000000000   0.0000000000
2     3.0000000000   0.0000000000   0.0000000000   1.0000000000
"""


def test_read_tmalign_matrix_handles_zero_and_one_based_rows(tmp_path: Path) -> None:
    zero_based = tmp_path / "m0.txt"
    zero_based.write_text(MATRIX_IDENTITY_SHIFT)
    one_based = tmp_path / "m1.txt"
    one_based.write_text(
        " ------ The rotation matrix to rotate Chain_1 to Chain_2 ------\n"
        "m               t[m]        u[m][0]        u[m][1]        u[m][2]\n"
        "1     1.0000000000   1.0000000000   0.0000000000   0.0000000000\n"
        "2     2.0000000000   0.0000000000   1.0000000000   0.0000000000\n"
        "3     3.0000000000   0.0000000000   0.0000000000   1.0000000000\n"
    )

    for matrix_path in (zero_based, one_based):
        translation, rotation = StructureOps.read_tmalign_matrix(matrix_path)
        assert translation.tolist() == [1.0, 2.0, 3.0]
        assert rotation.tolist() == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def test_transform_pdb_with_matrix_keeps_numbering_and_moves_coordinates(tmp_path: Path) -> None:
    """The whole point of the -m route: residue numbers survive the superposition."""
    source = tmp_path / "in.pdb"
    source.write_text(PDB_MULTI_CHAIN)
    matrix = tmp_path / "m.txt"
    matrix.write_text(MATRIX_IDENTITY_SHIFT)

    out = StructureOps.transform_pdb_with_matrix(source, tmp_path / "out.pdb", matrix)

    assert StructureOps.chain_ranges_from_pdb(out) == {"A": (5, 10), "B": (200, 210)}
    first_ca = next(
        line for line in out.read_text().splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CA"
    )
    x, y, z = float(first_ca[30:38]), float(first_ca[38:46]), float(first_ca[46:54])
    assert (x, y, z) == (pytest.approx(13.2), pytest.approx(14.3), pytest.approx(5.3))


def test_run_tmalign_matrix_superposes_structure_one(monkeypatch, tmp_path: Path) -> None:
    pdb1 = tmp_path / "a.pdb"
    pdb2 = tmp_path / "b.pdb"
    pdb1.write_text(PDB_MINI)
    pdb2.write_text(PDB_MINI)

    def _mock_run(cmd, cwd, capture_output, text, check):  # noqa: ANN001
        assert "-m" in cmd
        Path(cmd[cmd.index("-m") + 1]).write_text(MATRIX_IDENTITY_SHIFT)

        class _Result:
            stdout = "Aligned length= 2\n"

        return _Result()

    monkeypatch.setattr("common.structure_viz.subprocess.run", _mock_run)

    output = StructureOps.run_tmalign_matrix(pdb1, pdb2, tmp_path)
    assert output.reference == pdb2
    assert output.report == "Aligned length= 2\n"
    assert output.superposed.exists()
    # Structure 1 moved by the matrix's translation; the reference is untouched.
    assert StructureOps.chain_range_from_pdb(output.superposed, "A") == (1, 2)


def test_run_tmalign_matrix_raises_when_no_matrix_is_written(monkeypatch, tmp_path: Path) -> None:
    pdb1 = tmp_path / "a.pdb"
    pdb2 = tmp_path / "b.pdb"
    pdb1.write_text(PDB_MINI)
    pdb2.write_text(PDB_MINI)

    class _Result:
        stdout = ""

    monkeypatch.setattr(
        "common.structure_viz.subprocess.run",
        lambda cmd, cwd, capture_output, text, check: _Result(),
    )

    with pytest.raises(RuntimeError, match="no rotation matrix"):
        StructureOps.run_tmalign_matrix(pdb1, pdb2, tmp_path)


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


def test_tm_source_signature_prefers_upload_then_normalizes_pdb_id() -> None:
    upload = [{"datapath": "/tmp/file.pdb", "name": "2IVT.pdb"}]
    assert _tm_source_signature(upload, "1CRN") == ("upload", "/tmp/file.pdb|2IVT.pdb")
    assert _tm_source_signature(None, "2ivt") == ("pdb", "2IVT")
    assert _tm_source_signature(None, "") is None


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

    monkeypatch.setattr("common.structure_viz.requests.get", _mock_get)

    downloaded = service.download_pdb("2ivt")
    assert downloaded.name == "2IVT.pdb"
    assert downloaded.read_text() == PDB_MINI
