#!/usr/bin/env python3
"""
This script plots a QQ plot of observed vs expected -log10(p-values) for each cell type.

analysis-runner --dataset "bioheart" --description "plot qq plot" --access-level "test" \
    --output-dir "str/associatr/tob_n1055_and_bioheart_n990" --memory=8G \
    qqplotter.py \
    --input-dir=gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/raw_pval_extractor \
    --cell-types=CD4_TCM,CD4_Naive,NK,CD8_TEM,B_naive,CD8_Naive,CD14_Mono,CD4_TEM,CD8_TCM,B_intermediate,B_memory,Treg,CD4_CTL,gdT,CD16_Mono,MAIT,NK_CD56bright,cDC2,NK_Proliferating,dnT,pDC,Plasmablast,ILC,HSPC,CD8_Proliferating,cDC1,CD4_Proliferating,ASDC \
    --title='associaTR BioHEART' --ylim=315


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
    'CD4_TCM': '#0C46A0FF',
    'CD4_Naive': '#1976D2FF',
    'CD4_TEM': '#2096F2FF',
    'CD4_CTL': '#64B4F6FF',
    'Treg': '#90CAF8FF',
    'CD4_Proliferating': '#BADEFAFF',
    'gdT': '#817717FF',
    'MAIT': '#AEB32BFF',
    'dnT': '#CCDC39FF',
    'ILC': '#DCE674FF',
    'CD8_TEM': '#311A92FF',
    'CD8_Naive': '#5E34B1FF',
    'CD8_TCM': '#7E57C1FF',
    'CD8_Proliferating': '#D1C4E9FF',
    'NK': '#AC1357FF',
    'NK_CD56bright': '#E91E63FF',
    'NK_Proliferating': '#F38EB1FF',
    'B_naive': '#F47F17FF',
    'B_intermediate': '#FABF2CFF',
    'B_memory': '#FFEB3AFF',
    'Plasmablast': '#FFF176FF',
    'CD14_Mono': '#388D3BFF',
    'CD16_Mono': '#80C684FF',
    'cDC2': '#5D3F37FF',
    'pDC': '#795447FF',
    'cDC1': '#A0877FFF',
    'ASDC': '#D7CCC7FF',
    'HSPC': '#BDBDBDFF',
}

@click.option('--title', help='Title of the QQ plot')
@click.option('--ylim', help='Y-axis limit for the QQ plot', default=335)
@click.option('--input-dir', help='GCS path directory to the input gene-level p-value files')
@click.option('--cell-types', help='Comma-separated list of cell types to plot')
@click.command()
def main(input_dir, cell_types, title, ylim):
    init_batch()
    cell_type_list = cell_types.split(',')

    for cell_type in cell_type_list:
        df = pd.read_csv(
            f'{input_dir}/{cell_type}_gene_tests_raw_pvals.txt',
            header=None,
            names=['CHR', 'BP', 'raw_pval'],
            sep='\t',
        )
        df = df.dropna()

        globals()[f'observed_log_pvals_{cell_type}'] = -np.log10(df['raw_pval'])
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
    ax.set_ylabel('Expected -log₁₀(p-value)', fontsize=17)

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