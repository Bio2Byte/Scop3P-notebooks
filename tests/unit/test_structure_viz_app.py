from __future__ import annotations

import pandas as pd

from structure_viz.app import _b2b_table_dataframe, _selected_b2b_metric_column


def test_selected_b2b_metric_column_uses_toggle_state() -> None:
    assert _selected_b2b_metric_column("backbone", normalized=False) == "backbone"
    assert _selected_b2b_metric_column("backbone", normalized=True) == "backbone_normalized"
    assert _selected_b2b_metric_column(None, normalized=True) is None


def test_b2b_table_dataframe_switches_between_raw_and_normalized_values() -> None:
    dataframe = pd.DataFrame(
        {
            "Position": [1, 2],
            "Amino acid": ["A", "C"],
            "backbone": [0.1, 0.2],
            "sidechain": [0.3, 0.4],
            "ppII": [0.5, 0.6],
            "coil": [0.7, 0.8],
            "sheet": [0.9, 1.0],
            "helix": [1.1, 1.2],
            "earlyFolding": [1.3, 1.4],
            "disoMine": [1.5, 1.6],
            "backbone_normalized": [0.0, 1.0],
            "sidechain_normalized": [0.0, 1.0],
            "ppII_normalized": [0.0, 1.0],
            "coil_normalized": [0.0, 1.0],
            "sheet_normalized": [0.0, 1.0],
            "helix_normalized": [0.0, 1.0],
            "earlyFolding_normalized": [0.0, 1.0],
            "disoMine_normalized": [0.0, 1.0],
        }
    )

    raw_table = _b2b_table_dataframe(dataframe, normalized=False)
    normalized_table = _b2b_table_dataframe(dataframe, normalized=True)

    assert list(raw_table.columns) == [
        "Position",
        "Amino acid",
        "backbone",
        "sidechain",
        "ppII",
        "coil",
        "sheet",
        "helix",
        "earlyFolding",
        "disoMine",
    ]
    assert list(normalized_table.columns) == list(raw_table.columns)
    assert raw_table.loc[0, "backbone"] == 0.1
    assert normalized_table.loc[0, "backbone"] == 0.0
    assert raw_table.loc[1, "disoMine"] == 1.6
    assert normalized_table.loc[1, "disoMine"] == 1.0
