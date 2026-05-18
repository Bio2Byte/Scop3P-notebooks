from __future__ import annotations

import pandas as pd
import pytest

from common.mutation_effect import (
    Mutation,
    MutationEffectInference,
    MutationEffectService,
    MutationEffectViews,
)


def test_parse_mutations_accepts_csv_pairs() -> None:
    mutations = MutationEffectService.parse_mutations("10,25", "A,V")
    assert mutations == [Mutation(position=10, to_amino_acid="A"), Mutation(position=25, to_amino_acid="V")]


def test_apply_mutations_rewrites_sequence() -> None:
    sequence = "ACDEFG"
    mutated = MutationEffectService.apply_mutations(
        sequence,
        [Mutation(position=2, to_amino_acid="Y"), Mutation(position=6, to_amino_acid="W")],
    )
    assert mutated == "AYDEFW"


def test_build_mutation_labels_uses_wt_sequence() -> None:
    labels = MutationEffectService.build_mutation_labels(
        "ACDEFG",
        [Mutation(position=2, to_amino_acid="Y"), Mutation(position=6, to_amino_acid="W")],
    )
    assert labels == ["C2Y", "G6W"]


def test_inference_table_reports_label_shift() -> None:
    wt_df = pd.DataFrame(
        {
            "seqpos": [1, 2, 3],
            "seq": list("ABC"),
            "backbone": [0.70, 0.70, 0.70],
            "disoMine": [0.1, 0.2, 0.3],
            "earlyFolding": [0.1, 0.1, 0.1],
        }
    )
    mut_df = pd.DataFrame(
        {
            "seqpos": [1, 2, 3],
            "seq": list("AYC"),
            "backbone": [0.70, 1.10, 0.70],
            "disoMine": [0.1, 0.2, 0.3],
            "earlyFolding": [0.1, 0.1, 0.1],
        }
    )

    result = MutationEffectInference.mutation_effect_table_with_label_shift(
        wt_df=wt_df,
        mut_df=mut_df,
        feature="backbone",
        mutations=[Mutation(position=2, to_amino_acid="Y")],
        window=0,
    )

    assert result.iloc[0]["mutation"] == "B2Y"
    assert result.iloc[0]["label_shift_pos"] == "context-dependent -> membrane-spanning"


def test_make_wt_mut_merged_table_adds_ptm_flag() -> None:
    wt_df = pd.DataFrame(
        {
            "seqpos": [1, 2],
            "seq": list("AC"),
            "backbone": [0.7, 0.8],
            "disoMine": [0.2, 0.3],
            "earlyFolding": [0.1, 0.2],
        }
    )
    mut_df = pd.DataFrame(
        {
            "seqpos": [1, 2],
            "seq": list("AY"),
            "backbone": [0.7, 0.9],
            "disoMine": [0.2, 0.3],
            "earlyFolding": [0.1, 0.2],
        }
    )
    mods_df = pd.DataFrame({"position": [2]})

    merged = MutationEffectViews.make_wt_mut_merged_table(wt_df, mut_df, mods_df)
    assert list(merged.columns[:3]) == ["seqpos", "WT_AA", "Mut_AA"]
    assert merged.loc[merged["seqpos"] == 2, "PTMs"].iloc[0] == "yes"


class _TextResponse:
    def __init__(self, text: str) -> None:
        self.text = text

    def raise_for_status(self) -> None:
        return None


def test_fetch_uniprot_sequence_parses_fasta(monkeypatch) -> None:
    monkeypatch.setattr(
        "apps.common.mutation_effect.requests.get",
        lambda *args, **kwargs: _TextResponse(">sp|P12345|\nACD\nEFG\n"),
    )

    sequence = MutationEffectService().fetch_uniprot_sequence("P12345")
    assert sequence == "ACDEFG"


def test_fetch_uniprot_sequence_rejects_empty_sequence(monkeypatch) -> None:
    monkeypatch.setattr(
        "apps.common.mutation_effect.requests.get",
        lambda *args, **kwargs: _TextResponse(">sp|P12345|\n"),
    )

    with pytest.raises(ValueError, match="No sequence returned"):
        MutationEffectService().fetch_uniprot_sequence("P12345")


def test_fetch_scop3p_modifications_drops_invalid_positions(monkeypatch) -> None:
    class _Client:
        def fetch_modifications(self, accession: str) -> pd.DataFrame:
            return pd.DataFrame(
                [
                    {"position": 10, "residue": "S", "name": "Phosphorylation"},
                    {"position": pd.NA, "residue": "T", "name": "Phosphorylation"},
                    {"position": 12, "residue": "Y", "name": "Phosphorylation", "specificSinglyPhosphorylated": 1},
                ]
            )

    service = MutationEffectService(
        scop3p_client=_Client(),  # type: ignore[arg-type]
    )

    dataframe = service.fetch_scop3p_modifications("P12345")
    assert dataframe["position"].tolist() == [10, 12]
    assert list(dataframe.columns) == ["position", "residue", "name"]


def test_prediction_to_df_coerces_values_and_seqpos() -> None:
    prediction = {
        "proteins": {
            "P12345": {
                "seq": "AC",
                "backbone": ["0.7", "bad"],
                "disoMine": ["0.1", "0.2"],
                "earlyFolding": ["0.3", None],
            }
        }
    }

    dataframe = MutationEffectService.prediction_to_df(prediction, "P12345")
    assert dataframe["seqpos"].tolist() == [1, 2]
    assert dataframe.loc[0, "backbone"] == 0.7
    assert pd.isna(dataframe.loc[1, "backbone"])
    assert dataframe.loc[1, "disoMine"] == 0.2
