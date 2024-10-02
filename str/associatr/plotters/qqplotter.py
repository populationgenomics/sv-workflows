#!/usr/bin/env python3
"""
This script plots a QQ plot of observed vs expected -log10(p-values) for each cell type.

analysis-runner --dataset "bioheart" --description "plot qq plot" --access-level "test" \
    --output-dir "str/associatr/tob_n1055_and_bioheart_n990" --memory=8G \
    qqplotter.py \
    --input-dir=gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/raw_pval_extractor \
    --cell-types=CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,CD4_TCM_permuted,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC \
    --title='associaTR BioHEART' --ylim=330


"""
import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch, output_path


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

    # Define a list of colors
    colors = [
        'tab:blue',
        'tab:orange',
        'tab:green',
        'tab:red',
        'tab:purple',
        'tab:brown',
        'tab:pink',
        'tab:gray',
        'tab:olive',
        'tab:cyan',
        'darkred',
        'darkblue',
        'darkgreen',
        'darkorange',
        'darkviolet',
        'darkgoldenrod',
        'deeppink',
        'dimgray',
        'dodgerblue',
        'firebrick',
        'forestgreen',
        'fuchsia',
        'gold',
        'indigo',
        'khaki',
        'lawngreen',
        'lightcoral',
        'lightseagreen',
        'lightskyblue',
    ]

    # Pre-calculate sorted values
    expected_sorted_values = {
        cell_type: np.sort(globals()[f'expected_log_pvals_{cell_type}']) for cell_type in cell_type_list
    }
    observed_sorted_values = {
        cell_type: np.sort(globals()[f'observed_log_pvals_{cell_type}']) for cell_type in cell_type_list
    }

    # Loop through each cell type and plot the scatter plot
    for i, cell_type in enumerate(cell_type_list):
        output_label = cell_type_mapping.get(cell_type, cell_type)
        color_index = i % len(colors)  # Get the index of the color to use for the current cell type
        ax.scatter(
            expected_sorted_values[cell_type],
            observed_sorted_values[cell_type],
            color=colors[color_index],
            label=output_label,
            s=9,
        )

    # Create a legend for all items except "permuted control"
    handles, labels = ax.get_legend_handles_labels()
    permuted_control_handle = [h for h, l in zip(handles, labels) if l == "Permuted control"]
    other_handles = [h for h, l in zip(handles, labels) if l != "Permuted control"]
    other_labels = [l for l in labels if l != "Permuted control"]
    ax.add_artist(
        ax.legend(
            permuted_control_handle,
            ['Permuted control'],
            bbox_to_anchor=(1.05, 0),
            loc='upper left',
            fontsize=11,
        ),
    )
    sns.despine()
    # Create the main legend with other items
    ax.legend(other_handles, other_labels, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)

    ax.set_xlabel('Expected -log₁₀(p-value)', fontsize=14)
    ax.set_ylabel('Expected -log₁₀(p-value)', fontsize=14)

    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    ax.set_ylim(0, ylim)

    ax.plot([0, 7], [0, 7], color='grey', linestyle='--')  # Add a reference line

    gcs_output_path = output_path('summary_plots/publish/v1/qq_plot.png', 'analysis')
    fig.tight_layout()
    fig.savefig('qqplot.png')
    hl.hadoop_copy('qqplot.png', gcs_output_path)


if __name__ == '__main__':
    main()
