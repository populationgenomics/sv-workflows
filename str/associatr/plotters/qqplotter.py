#!/usr/bin/env python3
"""
This script plots a QQ plot of observed vs expected -log10(p-values) for each cell type.

analysis-runner --dataset "tenk10k" --description "plot qq plot" --access-level "test" \
    --output-dir "str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed/v2" --memory=32G \
    qqplotter.py \
    --input-dir=gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed/v2/raw_pval_extractor \
    --cell-types=CD4_TCM,CD4_Naive,NK,CD8_TEM,B_naive,CD14_Mono,CD8_Naive,CD4_TEM,CD8_TCM,B_intermediate,CD4_CTL,Treg,B_memory,CD16_Mono,MAIT,gdT,NK_CD56bright,cDC2,NK_Proliferating,dnT,pDC,Plasmablast,HSPC,ILC,CD4_Proliferating,CD8_Proliferating,cDC1,ASDC,CD4_TCM_permuted \
    --ylim=315
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import click
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import hail as hl

from cpg_utils.hail_batch import init_batch, output_path


# ---- Matplotlib style: prefer Arial, but fall back safely on HPC/containers ----
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


# Define the color mapping
color_mapping: Dict[str, str] = {
    "CD4_TCM": "#0C46A0FF",
    "CD4_Naive": "#1976D2FF",
    "CD4_TEM": "#2096F2FF",
    "CD4_CTL": "#64B4F6FF",
    "Treg": "#90CAF8FF",
    "CD4_Proliferating": "#BADEFAFF",
    "gdT": "#817717FF",
    "MAIT": "#AEB32BFF",
    "dnT": "#CCDC39FF",
    "ILC": "#DCE674FF",
    "CD8_TEM": "#311A92FF",
    "CD8_Naive": "#5E34B1FF",
    "CD8_TCM": "#7E57C1FF",
    "CD8_Proliferating": "#D1C4E9FF",
    "NK": "#AC1357FF",
    "NK_CD56bright": "#E91E63FF",
    "NK_Proliferating": "#F38EB1FF",
    "B_naive": "#F47F17FF",
    "B_intermediate": "#FABF2CFF",
    "B_memory": "#FFEB3AFF",
    "Plasmablast": "#FFF176FF",
    "CD14_Mono": "#388D3BFF",
    "CD16_Mono": "#80C684FF",
    "cDC2": "#5D3F37FF",
    "pDC": "#795447FF",
    "cDC1": "#A0877FFF",
    "ASDC": "#D7CCC7FF",
    "HSPC": "#BDBDBDFF",
}


# Cell type labels using mathtext for subscripts in legend
cell_type_mapping: Dict[str, str] = {
    "ASDC": "ASDC",
    "cDC1": "cDC1",
    "cDC2": "cDC2",
    "pDC": "pDC",
    "HSPC": "HSPC",
    "Plasmablast": "Plasmablast",
    "gdT": "gdT",
    "MAIT": "MAIT",
    "dnT": "dnT",
    "ILC": "ILC",
    "NK": "NK",
    "Treg": "Treg",
    # Subscripted labels
    "B_naive": r"$B_{\mathrm{Naive}}$",
    "B_memory": r"$B_{\mathrm{Memory}}$",
    "B_intermediate": r"$B_{\mathrm{Intermediate}}$",
    "NK_CD56bright": r"$NK_{\mathrm{CD56bright}}$",
    "NK_Proliferating": r"$NK_{\mathrm{Proliferating}}$",
    "CD14_Mono": r"$CD14_{\mathrm{Mono}}$",
    "CD16_Mono": r"$CD16_{\mathrm{Mono}}$",
    "CD4_Naive": r"$CD4_{\mathrm{Naive}}$",
    "CD4_TCM": r"$CD4_{\mathrm{TCM}}$",
    "CD4_TEM": r"$CD4_{\mathrm{TEM}}$",
    "CD4_CTL": r"$CD4_{\mathrm{CTL}}$",
    "CD4_Proliferating": r"$CD4_{\mathrm{Proliferating}}$",
    "CD8_Naive": r"$CD8_{\mathrm{Naive}}$",
    "CD8_TCM": r"$CD8_{\mathrm{TCM}}$",
    "CD8_TEM": r"$CD8_{\mathrm{TEM}}$",
    "CD8_Proliferating": r"$CD8_{\mathrm{Proliferating}}$",
    "CD4_TCM_permuted": "Permuted control",
}


def _read_pvals(input_dir: str, cell_type: str) -> pd.Series:
    """
    Reads raw p-values for a given cell type from a TSV in input_dir.

    Expects files named: <cell_type>_gene_tests_raw_pvals.txt
    with columns: CHR, BP, raw_pval (no header).
    """
    path = f"{input_dir}/{cell_type}_gene_tests_raw_pvals.txt"
    df = pd.read_csv(
        path,
        header=None,
        names=["CHR", "BP", "raw_pval"],
        sep="\t",
    ).dropna()

    # Guard against zeros (would become inf on -log10); also guard against >1.
    p = pd.to_numeric(df["raw_pval"], errors="coerce").dropna()
    p = p[(p > 0) & (p <= 1)]
    return p


def _expected_observed_logp(pvals: pd.Series) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns (expected_sorted, observed_sorted) for QQ plot in -log10 scale.
    """
    observed = -np.log10(pvals.to_numpy())
    n = observed.size
    expected = -np.log10(np.arange(1, n + 1) / n)
    return np.sort(expected), np.sort(observed)


@click.command()
@click.option("--title", help="Title of the QQ plot", default=None)
@click.option("--ylim", help="Y-axis limit for the QQ plot", default=335, type=float)
@click.option("--input-dir", help="GCS path directory to the input gene-level p-value files", required=True)
@click.option("--cell-types", help="Comma-separated list of cell types to plot", required=True)
def main(input_dir: str, cell_types: str, title: Optional[str], ylim: float) -> None:
    init_batch()

    cell_type_list: List[str] = [ct.strip() for ct in cell_types.split(",") if ct.strip()]
    if not cell_type_list:
        raise click.ClickException("No cell types provided via --cell-types")

    # Pre-calculate sorted values
    expected_sorted_values: Dict[str, np.ndarray] = {}
    observed_sorted_values: Dict[str, np.ndarray] = {}

    for cell_type in cell_type_list:
        pvals = _read_pvals(input_dir, cell_type)
        if pvals.empty:
            raise click.ClickException(f"No valid p-values found for cell type: {cell_type}")

        expected_sorted, observed_sorted = _expected_observed_logp(pvals)
        expected_sorted_values[cell_type] = expected_sorted
        observed_sorted_values[cell_type] = observed_sorted

    # Create QQ plot
    fig, ax = plt.subplots(figsize=(10, 8))

    default_color = "#808080"  # grey for unmapped cell types (incl permuted control)
    for cell_type in cell_type_list:
        output_label = cell_type_mapping.get(cell_type, cell_type)
        color = color_mapping.get(cell_type, default_color)

        ax.scatter(
            expected_sorted_values[cell_type],
            observed_sorted_values[cell_type],
            color=color,
            label=output_label,
            s=9,
            linewidths=0,
        )

    # Reference line
    max_x = max(float(vals.max()) for vals in expected_sorted_values.values())
    line_max = max(7.0, max_x)
    ax.plot([0, line_max], [0, line_max], color="grey", linestyle="--", linewidth=1)

    # Legends: split permuted control out
    handles, labels = ax.get_legend_handles_labels()
    permuted_control_idx = [i for i, l in enumerate(labels) if l == "Permuted control"]
    other_idx = [i for i, l in enumerate(labels) if l != "Permuted control"]

    if permuted_control_idx:
        ax.add_artist(
            ax.legend(
                [handles[permuted_control_idx[0]]],
                ["Permuted control"],
                bbox_to_anchor=(1.05, 0),
                loc="upper left",
                fontsize=12,
                frameon=False,
            )
        )

    other_handles = [handles[i] for i in other_idx]
    other_labels = [labels[i] for i in other_idx]
    ax.legend(
        other_handles,
        other_labels,
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        fontsize=11,
        frameon=False,
    )

    sns.despine()
    ax.set_xlabel("Expected -log₁₀(p-value)", fontsize=17)
    ax.set_ylabel("Observed -log₁₀(p-value)", fontsize=17)
    if title:
        ax.set_title(title, fontsize=17)

    ax.set_ylim(0, ylim)
    ax.tick_params(axis="both", labelsize=15)

    fig.tight_layout()

    # Save locally then copy to GCS
    local_out = "qqplot.png"
    fig.savefig(local_out, dpi=300)

    gcs_output_path = output_path("summary_plots/publish/v1/qq_plot.png", "analysis")
    hl.hadoop_copy(local_out, gcs_output_path)


if __name__ == "__main__":
    main()
