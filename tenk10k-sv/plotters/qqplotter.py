#!/usr/bin/env python3
"""
This script plots a QQ plot of observed vs expected -log10(p-values) for each cell type.

analysis-runner --dataset "tenk10k" --description "plot qq plot" --access-level "test" \
    --output-dir "tenk10k-sv/sv/bioheart_n968/7pc/v1" --memory=32G \
    --image 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:3c1041139794b7e91da004589cbbf2f592556faf-hail-95632ec3cb35a1e088b1604657e3eac8da853613' \
    qqplotter.py \
    --input-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/bioheart_n968/7pc/v1/raw_pval_extractor \
    --cell-types=CD4_TCM,CD4_TCM_permuted \
    --ylim=200


"""
import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch, output_path

# Define the color mapping
color_mapping = {
    'CD4_TCM': '#311A92FF',
    'CD4_Naive': '#512CA7FF',
    'CD4_TEM': '#6639B7FF',
    'CD4_CTL': '#9474CCFF',
    'Treg': '#B29DDAFF',
    'CD4_Proliferating': '#D1C4E9FF',
    'gdT': '#870D4EFF',
    'MAIT': '#C1185AFF',
    'dnT': '#E91E63FF',
    'ILC': '#F06192FF',
    'CD8_TEM': '#0C46A0FF',
    'CD8_Naive': '#1E87E5FF',
    'CD8_TCM': '#41A5F4FF',
    'CD8_Proliferating': '#BADEFAFF',
    'NK': '#2D7D32FF',
    'NK_CD56bright': '#4CAE50FF',
    'NK_Proliferating': '#A5D6A6FF',
    'B_naive': '#F47F17FF',
    'B_intermediate': '#FABF2CFF',
    'B_memory': '#FFEB3AFF',
    'Plasmablast': '#FFF176FF',
    'CD14_Mono': '#D84314FF',
    'CD16_Mono': '#FFAB91FF',
    'cDC2': '#5D3F37FF',
    'pDC': '#795447FF',
    'cDC1': '#A0877FFF',
    'ASDC': '#D7CCC7FF',
    'HSPC': '#BDBDBDFF',
}

def read_raw_pvals(input_dir, cell_type):
    """Read one cell type's raw p-values from raw_pval_extractor_sv.py output, as floats.

    The file is a headed 7-column TSV (chrom, pos, end, svtype, varid, gene, pval), so it is read
    with its header and the p-value is located by NAME -- the same convention as the other SV
    scripts. Reading it positionally instead (header=None + names=[...]) puts the literal string
    'pval' into the column, which makes it object/str dtype and makes -np.log10 raise
    'TypeError: loop of ufunc does not support argument 0 of type str'.
    """
    path = f'{input_dir}/{cell_type}_gene_tests_raw_pvals.tsv'
    with to_path(path).open() as handle:
        df = pd.read_csv(handle, sep='\t')

    if 'pval' not in df.columns:
        raise ValueError(f'{path} has no "pval" column; found {list(df.columns)}')

    # coerce rather than astype(float): this reports how many values were unusable
    # instead of raising on the first one
    pvals = pd.to_numeric(df['pval'], errors='coerce')

    n_bad = int(pvals.isna().sum())
    if n_bad:
        print(f'{cell_type}: dropping {n_bad:,} non-numeric/missing p-values')
    # subset to the p-value column only -- a frame-wide dropna() would also discard rows whose
    # svtype/end/varid are 'NA' (motif that did not parse), even though their p-value is fine
    pvals = pvals.dropna()

    # -log10(0) is inf and silently stretches the y-axis; associaTR can underflow to exactly 0
    n_zero = int((pvals == 0).sum())
    if n_zero:
        print(f'{cell_type}: {n_zero:,} p-values underflowed to 0, clamping to the smallest float')
        pvals = pvals.clip(lower=np.nextafter(0, 1))

    if pvals.dtype.kind != 'f':
        raise TypeError(f'{path}: pval read as {pvals.dtype}, not float')

    print(f'{cell_type}: {len(pvals):,} p-values from {path}')
    return pvals


@click.option('--title', help='Title of the QQ plot')
@click.option('--ylim', help='Y-axis limit for the QQ plot', default=335)
@click.option('--input-dir', help='GCS path directory to the input gene-level p-value files')
@click.option('--cell-types', help='Comma-separated list of cell types to plot')
@click.command()
def main(input_dir, cell_types, title, ylim):
    init_batch()
    cell_type_list = cell_types.split(',')

    for cell_type in cell_type_list:
        pvals = read_raw_pvals(input_dir, cell_type)

        globals()[f'observed_log_pvals_{cell_type}'] = -np.log10(pvals)
        globals()[f'n_{cell_type}'] = len(globals()[f'observed_log_pvals_{cell_type}'])
        globals()[f'expected_log_pvals_{cell_type}'] = -np.log10(
            np.arange(1, globals()[f'n_{cell_type}'] + 1) / globals()[f'n_{cell_type}'],
        )

    cell_type_mapping = {
        'ASDC': 'ASDC',
        'B_intermediate': 'B intermediate',
        'B_memory': 'B memory',
        'B_naive': 'B naive',
        'CD14_Mono': 'CD14+ Monocyte',
        'CD16_Mono': 'CD16+ Monocyte',
        'CD4_CTL': 'CD4+ CTL',
        'CD4_Naive': 'CD4+ Naive',
        'CD4_Proliferating': 'CD4+ Proliferating',
        'CD4_TCM': 'CD4+ TCM',
        'CD4_TCM_permuted': 'Permuted control',
        'CD4_TEM': 'CD4+ TEM',
        'CD8_Naive': 'CD8+ Naive',
        'CD8_Proliferating': 'CD8+ Proliferating',
        'CD8_TCM': 'CD8+ TCM',
        'CD8_TEM': 'CD8+ TEM',
        'cDC1': 'cDC1',
        'cDC2': 'cDC2',
        'dnT': 'dnT',
        'gdT': 'gdT',
        'HSPC': 'HSPC',
        'ILC': 'ILC',
        'MAIT': 'MAIT',
        'NK': 'NK',
        'NK_CD56bright': 'NK CD56bright',
        'NK_Proliferating': 'NK Proliferating',
        'pDC': 'pDC',
        'Plasmablast': 'Plasmablast',
        'Treg': 'Treg',
    }

    # Create QQ plot
    plt.figure(figsize=(10, 8))
    fig, ax = plt.subplots(figsize=(10, 8))

    # Set default color for permuted control or any cell type not in color_mapping
    default_color = '#808080'  # grey color for unmapped cell types

    # Pre-calculate sorted values
    expected_sorted_values = {
        cell_type: np.sort(globals()[f'expected_log_pvals_{cell_type}'])
        for cell_type in cell_type_list
    }
    observed_sorted_values = {
        cell_type: np.sort(globals()[f'observed_log_pvals_{cell_type}'])
        for cell_type in cell_type_list
    }

    # Plot each cell type in the order specified by cell_type_list
    for cell_type in cell_type_list:
        output_label = cell_type_mapping.get(cell_type, cell_type)
        # Use color from mapping if available, otherwise use default color
        color = color_mapping.get(cell_type, default_color)

        ax.scatter(
            expected_sorted_values[cell_type],
            observed_sorted_values[cell_type],
            color=color,
            label=output_label,
            s=9,
        )

    # Create a legend for permuted control and other items separately
    handles, labels = ax.get_legend_handles_labels()
    permuted_control_idx = [i for i, l in enumerate(labels) if l == "Permuted control"]
    other_idx = [i for i, l in enumerate(labels) if l != "Permuted control"]

    # Add permuted control legend if it exists
    if permuted_control_idx:
        ax.add_artist(
            ax.legend(
                [handles[permuted_control_idx[0]]],
                ["Permuted control"],
                bbox_to_anchor=(1.05, 0),
                loc='upper left',
                fontsize=12,
            )
        )

    # Create the main legend with other items, maintaining the order from cell_type_list
    other_handles = [handles[i] for i in other_idx]
    other_labels = [labels[i] for i in other_idx]
    ax.legend(other_handles, other_labels, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)

    sns.despine()
    ax.set_xlabel('Expected -log₁₀(p-value)', fontsize=17)
    ax.set_ylabel('Observed -log₁₀(p-value)', fontsize=17)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    ax.set_ylim(0, ylim)

    ax.plot([0, 7], [0, 7], color='grey', linestyle='--')  # Add a reference line

    gcs_output_path = output_path('summary_plots/publish/v1/qq_plot.png', 'analysis')
    fig.tight_layout()
    fig.savefig('qqplot.png')
    hl.hadoop_copy('qqplot.png', gcs_output_path)

if __name__ == '__main__':
    main()