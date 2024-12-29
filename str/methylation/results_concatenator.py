#!/usr/bin/env python3
"""
This script concatenates the results from running associatr-methylation on a per-chromosome basis.
The output is a single CSV file containing per chromosome results

analysis-runner --dataset bioheart --access-level test --output-dir potato --description "Concatenate associatr-methylation results" python3 results_concatenator.py \
    --chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,21,22


"""

import click
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

def process_file(file):
    df = pd.read_csv(file, sep='\t')
    site_chr = str(file).split('/')[-1].split('_')[0]
    site_pos = str(file).split('/')[-1].split('_')[1]
    site = f'{site_chr}_{site_pos}'
    old_p = f'p_{site}'
    old_coeff = f'coeff_{site}'
    old_se = f'se_{site}'

    # Rename the columns to be more generic
    df = df.rename(columns={old_p: 'p', old_coeff: 'coeff', old_se: 'se'})
    df['methyl_site'] = site  # add a column for the Cpg site coordinate
    return df

def concatenator(input_dir, chrom):
    files = list(to_path(f'{input_dir}/chr{chrom}').glob('*.tsv'))
    master_df = pd.DataFrame()

    # Use ThreadPoolExecutor to process files in parallel
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(process_file, files))

    # Concatenate all DataFrames
    master_df = pd.concat(results, ignore_index=True)
    output_gcs = f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr-methylation/bioheart_n25/5kb/results/chr{chrom}_concatenated_results.tsv'
    master_df.to_csv(output_gcs, sep='\t', index=False, header=True)


@click.option(
    '--input-dir',
    help='GCS path to the directory containing the associatr-methylation results',
    default='gs://cpg-bioheart-test-analysis/tenk10k/str/associatr-methylation/bioheart_n25/5kb/results',
)
@click.option('--chromosomes', help='Comma-separated list of chromosomes to concatenate', default='22')
@click.command()
def main(input_dir, chromosomes):
    b = get_batch(name='Concatenate associatr-methylation results')
    for chrom in chromosomes.split(','):
        j = b.new_python_job(
            name=f'Concatenate associatr-methylation results for chr{chrom}',
        )
        j.cpu(1)
        j.storage('5G')

        j.call(
            concatenator,
            input_dir,
            chrom,
        )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
