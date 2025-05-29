#!/usr/bin/env python3
"""
This script concatenates the results from running associatr-atac on a cell type basis.
The output is a single CSV file containing per cell type results

analysis-runner --dataset bioheart --access-level test --output-dir str/associatr-atac/tob/input_files/10kb_estrs/v1-remove-rare-GTs --description "Concatenate associatr-atac results" python3 results_concatenator.py \
    --celltypes "CD14_Mono,ASDC,CD14_Mono_permuted"


"""

import click
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

def process_file(file):
    df = pd.read_csv(file, sep='\t')
    site_pos = str(file).split('/')[-1].split('_')[1]
    site = str(file).split('/')[-1].split('_')[0]
    old_p = f'p_{site}'
    old_coeff = f'coeff_{site}'
    old_se = f'se_{site}'

    # Rename the columns to be more generic
    df = df.rename(columns={old_p: 'p', old_coeff: 'coeff', old_se: 'se'})
    df['atac_site'] = site  # add a column for the atac site coordinate
    return df

def concatenator(input_dir, cell_type):
    files = list(to_path(f'{input_dir}/{cell_type}').glob('*.tsv'))
    master_df = pd.DataFrame()

    # Use ThreadPoolExecutor to process files in parallel
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(process_file, files))

    # Concatenate all DataFrames
    master_df = pd.concat(results, ignore_index=True)
    output_gcs = f'gs://cpg-bioheart-test-analysis/str/associatr-atac/tob/input_files/10kb_estrs/v1-remove-rare-GTs/concat_results/{cell_type}_concatenated_results.tsv'
    master_df.to_csv(output_gcs, sep='\t', index=False, header=True)




@click.option(
    '--input-dir',
    help='GCS path to the directory containing the associatr-methylation results',
    default='gs://cpg-bioheart-test-analysis/str/associatr-atac/tob/input_files/10kb_estrs/v1-remove-rare-GTs/results',
)
@click.option('--celltypes', help='Comma-separated list of cell types to concatenate')
@click.command()
def main(input_dir, celltypes):
    b = get_batch(name='Concatenate associatr-atac results')
    for cell_type in celltypes.split(','):
        j = b.new_python_job(
            name=f'Concatenate associatr-atac results for {cell_type}',
        )
        j.cpu(1)
        j.storage('5G')

        j.call(
            concatenator,
            input_dir,
            cell_type,
        )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter