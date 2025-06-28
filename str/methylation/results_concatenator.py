#!/usr/bin/env python3
"""
This script concatenates the results from running associatr-methylation on a per-chromosome basis.
The output is a single CSV file containing per chromosome results

analysis-runner --dataset tenk10k --access-level test --output-dir potato --description "Concatenate associatr-methylation results" python3 results_concatenator.py \
    --chromosomes 22


"""

import click
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

def process_file(file):
    try:
        df = pd.read_csv(file, sep='\t')
        site_chr = str(file).split('/')[-1].split('_')[0]
        site_pos = str(file).split('/')[-1].split('_')[1]
        site = f'{site_chr}_{site_pos}'
        old_p = f'p_{site}'
        old_coeff = f'coeff_{site}'
        old_se = f'se_{site}'

        # Rename the columns to be more generic
        df = df.rename(columns={old_p: 'p', old_coeff: 'coeff', old_se: 'se'})
        df['methyl_site'] = site
        return df
    except Exception as e:
        print(f"Failed processing {file}: {e}")
        return None

def concatenator(input_dir, chrom):
    files = list(to_path(f'{input_dir}/chr{chrom}').glob('*.tsv'))
    dataframes = []

    for file in files:
        df = process_file(file)
        if df is not None:
            dataframes.append(df)

    if dataframes:
        master_df = pd.concat(dataframes, ignore_index=True)
        output_gcs = f'gs://cpg-tenk10k-test-analysis/str/associatr-methylation/bioheart_n25/5kb/results/chr{chrom}_concatenated_results.tsv'
        master_df.to_csv(output_gcs, sep='\t', index=False, header=True)
        print(f"Wrote {len(master_df)} rows to {output_gcs}")
    else:
        print("No valid dataframes to concatenate.")




@click.option(
    '--input-dir',
    help='GCS path to the directory containing the associatr-methylation results',
    default='gs://cpg-tenk10k-test-analysis/str/associatr-methylation/bioheart_n25/5kb/results',
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
