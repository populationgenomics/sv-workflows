#!/usr/bin/env python3

"""
This script extracts the raw p-values from the results of associaTR into one text file per cell type.
For downstream use to make a QQ plot.

analysis-runner --dataset "bioheart" --description "raw pval extractor" --access-level "test" \
    --output-dir "str/associatr/240_libraries_tenk10kp1_v2_run/results" \
    raw_pval_extractor.py --input-dir=gs://cpg-bioheart-test/str/associatr/240_libraries_tenk10kp1_v2_run/results/v1 \
    --cell-types=CD8_TEM --chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22

"""

import pandas as pd

import click
from cpg_utils.hail_batch import output_path
from cpg_utils import to_path


@click.option(
    '--input-dir',
    help='GCS path to the raw results of associaTR',
    type=str,
)
@click.option(
    '--cell-types',
    help='Name of the cell type, comma separated if multiple',
)
@click.option(
    '--chromosomes',
    help='Chromosome number eg 1, comma separated if multiple',
)
@click.command()
def main(input_dir, cell_types, chromosomes):

    for cell_type in cell_types.split(','):
        gcs_output = output_path(
            f'raw_pval_extractor/{cell_type}_gene_tests_raw_pvals.txt'
        )
        with to_path(gcs_output).open('w') as f:
            for chromosome in chromosomes.split(','):
                gene_files = list(
                    to_path(f'{input_dir}/{cell_type}/chr{chromosome}').glob('*.tsv')
                )
                for gene_file in gene_files:
                    # read the raw results
                    gene_results = pd.read_csv(gene_file, sep='\t')
                    pvals = gene_results.iloc[:, 5]
                    for pval in pvals:
                        f.write(str(pval) + '\n')


if __name__ == '__main__':
    main()
