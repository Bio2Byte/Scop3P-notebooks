from __future__ import annotations

from pathlib import Path

from common.session_metadata import build_metadata, write_metadata


def test_build_metadata_contains_context_only_fields(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("SCOP3P_APP_NAME", "peptide-mapper")
    monkeypatch.setenv("SCOP3P_IMAGE_VERSION", "v9")
    monkeypatch.setenv("SCOP3P_IMAGE_REVISION", "deadbeef")

    log_file = tmp_path / "scop3p_toolkit_log_20260519_120000_000001.log"
    payload = build_metadata(
        log_file_path=log_file,
        log_dir=tmp_path,
        session_started_at="2026-05-19T12:00:00+00:00",
    )

    assert payload["application"]["name"] == "peptide-mapper"
    assert payload["session"]["started_at_utc"] == "2026-05-19T12:00:00+00:00"
    assert payload["session"]["log_file"] == str(log_file)
    assert payload["image"]["version"] == "v9"
    assert payload["image"]["revision"] == "deadbeef"
    assert "python" in payload["runtime"]
    assert "shiny" in payload["dependencies"]
    assert "TM-align" in payload["external_tools"]


def test_write_metadata_creates_yaml_file(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("SCOP3P_APP_NAME", "structure-viz")
    metadata_path = tmp_path / "metadata.yml"
    log_file = tmp_path / "scop3p_toolkit_log_20260519_120000_000001.log"

    result = write_metadata(
        metadata_path=metadata_path,
        log_file_path=log_file,
        log_dir=Path(tmp_path),
        session_started_at="2026-05-19T12:00:00+00:00",
    )

    assert result == metadata_path
    contents = metadata_path.read_text(encoding="utf-8")
    assert 'schema_version: "1.0"' in contents
    assert 'name: "structure-viz"' in contents
    assert f'log_file: "{log_file}"' in contents
