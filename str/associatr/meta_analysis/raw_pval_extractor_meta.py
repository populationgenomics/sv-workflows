#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
This script extracts the raw p-values from the results of the meta-analysis output files.
For downstream use to make a QQ plot.

analysis-runner --dataset "bioheart" --description "raw pval extractor" --access-level "test" \
    --output-dir "str/associatr/bioheart_n990/EUR_only/results" \
    raw_pval_extractor.py --input-dir=gs://cpg-bioheart-test-analysis/str/associatr/bioheart_n990/EUR_only/results/v1 \
    --cell-types=CD4_TCM_permuted --chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22

"""

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import output_path


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
    """
    Extracts the raw p-values from the results of associaTR into one text file per cell type.
    """
    for cell_type in cell_types.split(','):
        gcs_output = output_path(f'raw_pval_extractor/{cell_type}_gene_tests_raw_pvals.txt')
        with to_path(gcs_output).open('w') as f:
            for chromosome in chromosomes.split(','):
                gene_files = list(to_path(f'{input_dir}/{cell_type}/chr{chromosome}').glob('*.tsv'))
                for gene_file in gene_files:
                    # read the raw results
                    gene_results = pd.read_csv(gene_file, sep='\t')
                    chr = gene_results.iloc[:, 0]
                    pos = gene_results.iloc[:, 1]
                    pvals = gene_results.iloc[:, 4]
                    for chr1, pos1, pval1 in zip(chr, pos, pvals):
                        f.write(chr1 + '\t' + str(pos1) + '\t' + str(pval1) + '\n')


if __name__ == '__main__':
    main()
