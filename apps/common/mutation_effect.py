from __future__ import annotations

from dataclasses import dataclass
from tempfile import NamedTemporaryFile
from uuid import uuid4

import pandas as pd
import requests
from bokeh.embed import components
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from bokeh.resources import CDN

from common.services import Scop3PClient

try:
    from b2bTools import SingleSeq
    from b2bTools import constants
except Exception:  # pragma: no cover - exercised in container/runtime
    SingleSeq = None
    constants = None


VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(slots=True, frozen=True)
class Mutation:
    position: int
    to_amino_acid: str

    @property
    def label(self) -> str:
        return f"{self.position}{self.to_amino_acid}"


class MutationEffectService:
    """Shared logic for the mutation-effect notebook conversion."""

    def __init__(
        self,
        scop3p_client: Scop3PClient | None = None,
        uniprot_base_url: str = "https://rest.uniprot.org/uniprotkb",
        timeout: int = 30,
    ) -> None:
        self.scop3p_client = scop3p_client or Scop3PClient(timeout=timeout)
        self.uniprot_base_url = uniprot_base_url.rstrip("/")
        self.timeout = timeout

    def fetch_uniprot_sequence(self, accession: str) -> str:
        response = requests.get(f"{self.uniprot_base_url}/{accession}.fasta", timeout=self.timeout)
        response.raise_for_status()
        sequence = "".join(
            line.strip()
            for line in response.text.splitlines()
            if line and not line.startswith(">")
        )
        if not sequence:
            raise ValueError(f"No sequence returned for {accession}.")
        return sequence

    def fetch_scop3p_modifications(self, accession: str) -> pd.DataFrame:
        dataframe = self.scop3p_client.fetch_modifications(accession)
        if dataframe.empty:
            return dataframe

        keep = [column for column in dataframe.columns if column != "specificSinglyPhosphorylated"]
        dataframe = dataframe[keep].copy()
        dataframe["position"] = pd.to_numeric(dataframe["position"], errors="coerce")
        dataframe = dataframe.dropna(subset=["position"])
        dataframe["position"] = dataframe["position"].astype(int)
        return dataframe

    def predict_biophysical(self, accession: str, sequence: str) -> dict:
        if SingleSeq is None:
            raise RuntimeError("b2bTools is not available in this environment.")

        with NamedTemporaryFile(prefix="seq_", suffix=".fasta", mode="w") as handle:
            handle.write(f">{accession}\n{sequence}\n")
            handle.flush()
            predictor = SingleSeq(handle.name)
            
            predictor.predict(
                tools=[constants.TOOL_DYNAMINE, constants.TOOL_DISOMINE, constants.TOOL_EFOLDMINE]
            )
            
            return predictor.get_all_predictions()

    @staticmethod
    def prediction_to_df(prediction: dict, accession: str) -> pd.DataFrame:
        protein = None
        if isinstance(prediction, dict):
            if "proteins" in prediction and accession in prediction["proteins"]:
                protein = prediction["proteins"][accession]
            elif accession in prediction:
                protein = prediction[accession]
        if protein is None:
            raise ValueError("Could not find predictions for the requested accession.")

        size = (
            len(protein.get("seq", ""))
            or len(protein.get("backbone", []))
            or len(protein.get("disoMine", []))
            or len(protein.get("earlyFolding", []))
        )
        dataframe = pd.DataFrame(
            {
                "seq": protein.get("seq", [None] * size),
                "seqpos": list(range(1, size + 1)),
                "backbone": protein.get("backbone", [None] * size),
                "sidechain": protein.get("sidechain", [None] * size),
                "helix": protein.get("helix", [None] * size),
                "sheet": protein.get("sheet", [None] * size),
                "coil": protein.get("coil", [None] * size),
                "ppII": protein.get("ppII", [None] * size),
                "disoMine": protein.get("disoMine", [None] * size),
                "earlyFolding": protein.get("earlyFolding", [None] * size),
            }
        )
        for column in ["backbone", "sidechain", "helix", "sheet", "coil", "ppII", "disoMine", "earlyFolding"]:
            dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")

        return dataframe

    @staticmethod
    def parse_mutations(positions_csv: str, aas_csv: str) -> list[Mutation]:
        positions = [value.strip() for value in positions_csv.split(",") if value.strip()]
        amino_acids = [value.strip().upper() for value in aas_csv.split(",") if value.strip()]
        if len(positions) != len(amino_acids):
            raise ValueError("Number of positions must match number of amino acids.")

        mutations: list[Mutation] = []
        for position_raw, amino_acid in zip(positions, amino_acids, strict=True):
            if not position_raw.isdigit():
                raise ValueError(f"Invalid position '{position_raw}'.")
            position = int(position_raw)
            if len(amino_acid) != 1 or amino_acid not in VALID_AMINO_ACIDS:
                raise ValueError(f"Invalid amino acid '{amino_acid}' at position {position}.")
            mutations.append(Mutation(position=position, to_amino_acid=amino_acid))
        return mutations

    @staticmethod
    def apply_mutations(sequence: str, mutations: list[Mutation]) -> str:
        residues = list(sequence)
        for mutation in mutations:
            if mutation.position < 1 or mutation.position > len(residues):
                raise ValueError(
                    f"Position {mutation.position} is out of range for sequence length {len(residues)}."
                )
            residues[mutation.position - 1] = mutation.to_amino_acid
        return "".join(residues)

    @staticmethod
    def build_mutation_labels(sequence: str, mutations: list[Mutation]) -> list[str]:
        labels = []
        for mutation in mutations:
            if 1 <= mutation.position <= len(sequence):
                labels.append(f"{sequence[mutation.position - 1]}{mutation.position}{mutation.to_amino_acid}")
            else:
                labels.append(f"?{mutation.position}{mutation.to_amino_acid}")
        return labels


class MutationEffectInference:
    """Pure helpers for classification labels and summary tables."""

    @staticmethod
    def label_backbone(value: float) -> str:
        if pd.isna(value):
            return "NA"
        if value > 1.0:
            return "membrane-spanning"
        if value > 0.8:
            return "rigid"
        if value > 0.69:
            return "context-dependent"
        return "flexible"

    @staticmethod
    def label_disorder(value: float) -> str:
        if pd.isna(value):
            return "NA"
        return "disordered" if value > 0.50 else "ordered"

    @staticmethod
    def label_earlyfold(value: float) -> str:
        if pd.isna(value):
            return "NA"
        return "early-folding" if value > 0.169 else "non-early-folding"

    LABEL_FUNCS = {
        "backbone": (label_backbone.__func__, "Backbone dynamics"),
        "disoMine": (label_disorder.__func__, "Disorder"),
        "earlyFolding": (label_earlyfold.__func__, "Early folding"),
    }

    @classmethod
    def mutation_effect_table_with_label_shift(
        cls,
        wt_df: pd.DataFrame,
        mut_df: pd.DataFrame,
        feature: str,
        mutations: list[Mutation],
        window: int = 5,
    ) -> pd.DataFrame:
        if feature not in cls.LABEL_FUNCS:
            raise ValueError(f"feature must be one of {list(cls.LABEL_FUNCS.keys())}")

        label_fn, pretty = cls.LABEL_FUNCS[feature]
        wt = wt_df.set_index("seqpos", drop=False)
        mut = mut_df.set_index("seqpos", drop=False)
        max_pos = int(wt["seqpos"].max())

        rows = []
        for mutation in mutations:
            if mutation.position not in wt.index or mutation.position not in mut.index:
                continue

            wt_center = float(wt.loc[mutation.position, feature])
            mut_center = float(mut.loc[mutation.position, feature])
            delta_center = mut_center - wt_center

            low = max(1, mutation.position - window)
            high = min(max_pos, mutation.position + window)
            wt_mean = float(wt.loc[low:high, feature].astype(float).mean())
            mut_mean = float(mut.loc[low:high, feature].astype(float).mean())
            delta_mean = mut_mean - wt_mean

            wt_label_center = label_fn(wt_center)
            mut_label_center = label_fn(mut_center)
            wt_label_mean = label_fn(wt_mean)
            mut_label_mean = label_fn(mut_mean)

            shift_center = (
                f"{wt_label_center} -> {mut_label_center}"
                if wt_label_center != mut_label_center
                else f"{wt_label_center} (no class change)"
            )
            shift_mean = (
                f"{wt_label_mean} -> {mut_label_mean}"
                if wt_label_mean != mut_label_mean
                else f"{wt_label_mean} (no class change)"
            )

            wt_aa = str(wt.loc[mutation.position, "seq"])
            mut_aa_seq = str(mut.loc[mutation.position, "seq"])
            notes = []
            if mut_aa_seq.upper() != mutation.to_amino_acid.upper():
                notes.append(f"Mut seq AA={mut_aa_seq} (expected {mutation.to_amino_acid})")

            rows.append(
                {
                    "pos": mutation.position,
                    "WT_AA": wt_aa,
                    "Mut_AA": mutation.to_amino_acid,
                    "mutation": f"{wt_aa}{mutation.position}{mutation.to_amino_acid}",
                    f"{feature}_WT@pos": wt_center,
                    f"{feature}_Mut@pos": mut_center,
                    "delta_pos": delta_center,
                    f"{feature}_WT_mean": wt_mean,
                    f"{feature}_Mut_mean": mut_mean,
                    "delta_mean": delta_mean,
                    "label_shift_pos": shift_center,
                    "label_shift_mean": shift_mean,
                    "inference": (
                        f"{pretty}: {shift_center} at site (delta {delta_center:+.3f}); "
                        f"window mean: {shift_mean} (delta {delta_mean:+.3f})"
                    ),
                    "note": "; ".join(notes),
                }
            )

        dataframe = pd.DataFrame(rows)
        numeric_columns = [
            column
            for column in dataframe.columns
            if any(key in column for key in ["_WT@pos", "_Mut@pos", "_WT_mean", "_Mut_mean", "delta_"])
        ]
        for column in numeric_columns:
            dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce").round(3)
        return dataframe


class MutationEffectViews:
    """Rendering helpers for Bokeh plots, guides, and scrollable tables."""

    @staticmethod
    def bokeh_document(fig) -> str:  # noqa: ANN001
        script, div = components(fig)
        return CDN.render() + div + script

    @staticmethod
    def make_wt_plot(wt_df: pd.DataFrame, mods_df: pd.DataFrame) -> str:
        plot = figure(
            width=1000,
            height=300,
            tools="pan,box_zoom,reset,save",
            toolbar_location="below",
            toolbar_sticky=False,
        )
        plot.title.text = "Biophysical properties (WT)"

        backbone = plot.line(
            wt_df["seqpos"],
            wt_df["backbone"],
            line_width=2,
            color="blue",
            alpha=0.8,
            muted_color="blue",
            muted_alpha=0.2,
            legend_label="backbone_dynamics",
        )
        disorder = plot.line(
            wt_df["seqpos"],
            wt_df["disoMine"],
            line_width=2,
            color="red",
            alpha=0.8,
            muted_color="red",
            muted_alpha=0.2,
            legend_label="disorder",
        )
        early = plot.line(
            wt_df["seqpos"],
            wt_df["earlyFolding"],
            line_width=2,
            color="grey",
            alpha=0.8,
            muted_color="grey",
            muted_alpha=0.2,
            legend_label="earlyFolding",
        )
        plot.add_tools(HoverTool(tooltips="Seqpos:@x, value:@y", renderers=[backbone, disorder, early]))
        MutationEffectViews._add_ptm_markers(plot, mods_df)
        plot.legend.click_policy = "mute"
        plot.add_layout(plot.legend[0], "right")
        return MutationEffectViews.bokeh_document(plot)

    @staticmethod
    def make_mut_plot(wt_df: pd.DataFrame, mut_df: pd.DataFrame, mods_df: pd.DataFrame) -> str:
        plot = figure(
            width=1000,
            height=300,
            tools="pan,box_zoom,reset,save",
            toolbar_location="below",
            toolbar_sticky=False,
        )
        plot.title.text = "Biophysical properties (WT vs Mutant)"

        renderers = [
            plot.line(wt_df["seqpos"], wt_df["backbone"], line_width=2, color="skyblue", alpha=0.8, muted_color="skyblue", muted_alpha=0.2, legend_label="backbone_dynamics (WT)"),
            plot.line(wt_df["seqpos"], wt_df["disoMine"], line_width=2, color="salmon", alpha=0.8, muted_color="salmon", muted_alpha=0.2, legend_label="disorder (WT)"),
            plot.line(wt_df["seqpos"], wt_df["earlyFolding"], line_width=2, color="grey", alpha=0.8, muted_color="grey", muted_alpha=0.2, legend_label="earlyFolding (WT)"),
            plot.line(mut_df["seqpos"], mut_df["backbone"], line_width=2, color="blue", alpha=0.8, muted_color="blue", muted_alpha=0.2, legend_label="backbone_mut"),
            plot.line(mut_df["seqpos"], mut_df["disoMine"], line_width=2, color="red", alpha=0.8, muted_color="red", muted_alpha=0.2, legend_label="disorder_mut"),
            plot.line(mut_df["seqpos"], mut_df["earlyFolding"], line_width=2, color="black", alpha=0.8, muted_color="black", muted_alpha=0.2, legend_label="earlyFolding_mut"),
        ]
        plot.add_tools(HoverTool(tooltips="Seqpos:@x, value:@y", renderers=renderers))
        MutationEffectViews._add_ptm_markers(plot, mods_df)
        plot.legend.click_policy = "mute"
        plot.add_layout(plot.legend[0], "right")
        return MutationEffectViews.bokeh_document(plot)

    @staticmethod
    def _add_ptm_markers(plot, mods_df: pd.DataFrame) -> None:  # noqa: ANN001
        if mods_df is None or mods_df.empty or "position" not in mods_df.columns:
            return

        source = ColumnDataSource(
            dict(
                x=mods_df["position"].tolist(),
                y=[0.5] * len(mods_df),
                residue=mods_df.get("residue", pd.Series([""] * len(mods_df))).astype(str).tolist(),
                name=mods_df.get("name", pd.Series([""] * len(mods_df))).astype(str).tolist(),
                source=mods_df.get("source", pd.Series([""] * len(mods_df))).astype(str).tolist(),
            )
        )
        renderer = plot.scatter(
            x="x",
            y="y",
            source=source,
            marker="circle",
            size=10,
            fill_alpha=0.6,
            line_alpha=0.8,
            color="grey",
            legend_label="P-sites",
        )
        plot.add_tools(
            HoverTool(
                tooltips=[("Seqpos", "@x"), ("residue", "@residue"), ("mod", "@name"), ("source", "@source")],
                renderers=[renderer],
            )
        )

    @staticmethod
    def track_guide_html() -> str:
        return """
<div style="margin:6px 0 10px 0; padding:8px 10px; border:1px solid #ddd; border-radius:8px;">
  <div style="margin-bottom:6px;"><b>How to read the tracks</b></div>
  <div style="margin:2px 0;"><span style="color:blue; font-weight:600;">Backbone dynamics</span>: &gt;1.0 membrane-spanning, 0.8-1.0 rigid, 0.69-0.80 context-dependent, &lt;0.69 flexible</div>
  <div style="margin:2px 0;"><span style="color:red; font-weight:600;">Disorder (DisoMine)</span>: values &gt;0.50 indicate disordered regions</div>
  <div style="margin:2px 0;"><span style="color:grey; font-weight:600;">Early folding</span>: values &gt;0.169 suggest early-folding propensity</div>
  <div style="margin:2px 0;"><span style="color:grey; font-weight:600;">P-sites</span>: phosphorylation positions (grey dots)</div>
</div>
"""

    @staticmethod
    def make_wt_table(wt_df: pd.DataFrame, mods_df: pd.DataFrame) -> pd.DataFrame:
        columns = MutationEffectViews._non_runtime_pred_cols(wt_df)
        base = wt_df[["seqpos", "seq"] + columns].copy()
        return MutationEffectViews._add_ptm_flag(base, mods_df)

    @staticmethod
    def make_wt_mut_merged_table(
        wt_df: pd.DataFrame,
        mut_df: pd.DataFrame,
        mods_df: pd.DataFrame,
    ) -> pd.DataFrame:
        wt_columns = MutationEffectViews._non_runtime_pred_cols(wt_df)
        mut_columns = MutationEffectViews._non_runtime_pred_cols(mut_df)
        common = [column for column in wt_columns if column in mut_columns]
        wt = wt_df[["seqpos", "seq"] + common].copy().rename(columns={"seq": "WT_AA"})
        mut = mut_df[["seqpos", "seq"] + common].copy().rename(columns={"seq": "Mut_AA"})
        merged = pd.merge(wt, mut, on="seqpos", how="inner", suffixes=("_WT", "_Mut"))

        ordered = ["seqpos", "WT_AA", "Mut_AA"]
        for column in common:
            ordered.extend([f"{column}_WT", f"{column}_Mut"])
        merged = merged[ordered]
        return MutationEffectViews._add_ptm_flag(merged, mods_df)

    @staticmethod
    def _non_runtime_pred_cols(dataframe: pd.DataFrame) -> list[str]:
        core = {"seqpos", "seq"}
        return [
            column
            for column in dataframe.columns
            if column not in core and "runtime" not in str(column).lower()
        ]

    @staticmethod
    def _add_ptm_flag(dataframe: pd.DataFrame, mods_df: pd.DataFrame) -> pd.DataFrame:
        output = dataframe.copy()
        positions = set()
        if mods_df is not None and not mods_df.empty and "position" in mods_df.columns:
            positions = set(pd.to_numeric(mods_df["position"], errors="coerce").dropna().astype(int).tolist())
        output["PTMs"] = output["seqpos"].apply(lambda value: "yes" if int(value) in positions else "no")
        return output

    @staticmethod
    def scrollable_table_html(
        dataframe: pd.DataFrame,
        *,
        title: str = "",
        height_px: int = 420,
        width: str = "100%",
        sticky_cols: int = 0,
        col_widths_px: list[int] | None = None,
        highlight_seqpos: list[int] | None = None,
        seqpos_col: str = "seqpos",
    ) -> str:
        if dataframe is None or dataframe.empty:
            return f"<div><b>{title}</b><br><i>(no rows)</i></div>" if title else "<i>(no rows)</i>"

        table_id = f"tbl_{uuid4().hex}"
        if col_widths_px is None:
            col_widths_px = [90] * sticky_cols
        else:
            col_widths_px = (list(col_widths_px) + [90] * sticky_cols)[:sticky_cols]

        left_offsets = []
        offset = 0
        for width_px in col_widths_px:
            left_offsets.append(offset)
            offset += int(width_px)

        highlight = set(int(value) for value in highlight_seqpos) if highlight_seqpos else set()
        if highlight and seqpos_col in dataframe.columns:
            columns = list(dataframe.columns)
            thead = "<thead><tr>" + "".join(f"<th>{column}</th>" for column in columns) + "</tr></thead>"
            rows = []
            for _, row in dataframe.iterrows():
                try:
                    seqpos = int(row[seqpos_col])
                except Exception:
                    seqpos = None
                row_class = " class='row-hl'" if seqpos is not None and seqpos in highlight else ""
                tds = "".join(f"<td>{'' if pd.isna(row[column]) else str(row[column])}</td>" for column in columns)
                rows.append(f"<tr{row_class}>{tds}</tr>")
            html_table = f"<table border='1' class='dataframe'>{thead}<tbody>{''.join(rows)}</tbody></table>"
        else:
            html_table = dataframe.to_html(index=False, escape=True)

        sticky_css = []
        for index in range(1, sticky_cols + 1):
            sticky_css.append(
                f"""
#{table_id} table th {{ text-align:center;}},
#{table_id} table th:nth-child({index}),
#{table_id} table td:nth-child({index}) {{
  position: sticky;
  left: {left_offsets[index - 1]}px;
  z-index: {4 + (sticky_cols - index)};
  background: #fff;
}}
"""
            )

        colgroup = "<colgroup>" + "".join(
            f"<col style='width:{int(width_px)}px;'>"
            for width_px in col_widths_px
        ) + "</colgroup>"
        html_table = html_table.replace(
            "<table border=\"1\" class=\"dataframe\">",
            f"<table border=\"1\" class=\"dataframe\">{colgroup}",
            1,
        )

        return f"""
<style>
#{table_id} .tbl-title {{ font-weight: 700; margin: 8px 0 6px 0; }}
#{table_id} tr.row-hl td {{ background: #fff6cc; }}
#{table_id} .tbl-wrap {{
  width: {width};
  max-width: {width};
  height: {height_px}px;
  overflow: auto;
  border: 1px solid #ddd;
  border-radius: 8px;
}}
#{table_id} table {{ border-collapse: collapse; width: max-content; min-width: 100%; font-size: 12px; }}
#{table_id} th, #{table_id} td {{ border: 1px solid #eee; padding: 6px 8px; white-space: nowrap; }}
#{table_id} thead th {{ position: sticky; top: 0; z-index: 10; background: #f7f7f7; }}
{''.join(sticky_css)}
</style>
<div id="{table_id}">
  {f'<div class="tbl-title">{title}</div>' if title else ''}
  <div class="tbl-wrap">{html_table}</div>
</div>
"""
