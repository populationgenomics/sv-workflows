#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""

This script performs FDR (across gene) correction (the second and final step of multiple testing correction).
Ensure that `run_acat.py` has been run to generate the gene-level p-values first.
Output is one TSV file per cell type with three columns - gene name, gene-level p-value (ACAT-corrected), and q-value.

analysis-runner --dataset "bioheart" --description "compute qvals" --access-level "test" \
    --output-dir "str/associatr/rna_pc_calibration/2_pcs/results" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images-dev/r-qvalue:1.0 \
    run_storey.py --input-dir=gs://cpg-bioheart-test/str/associatr/rna_pc_calibration/2_pcs/results/gene_level_pvals \
    --cell-types=CD8_TEM --chromosomes=2

"""

import pandas as pd

import click
from cpg_utils.hail_batch import get_batch, output_path
from cpg_utils import to_path
import rpy2.robjects as ro


def compute_storey(input_dir, cell_type, chromosomes):
    """
    Compute Storey's q-values for gene-level p-values
    """
    ro.r('library(qvalue)')

    for chromosome in chromosomes.split(','):
        # read in gene-level p-values
        gene_pval_files = list(
            to_path(f'{input_dir}/{cell_type}/chr{chromosome}').glob('*.tsv')
        )
        pval_df = pd.read_csv(gene_pval_files[0], sep='\t')
        for gene_pval_file in gene_pval_files[1:]:
            pval_df = pd.concat([pval_df, pd.read_csv(gene_pval_file, sep='\t')])

    pvals = pval_df['gene_level_pval']
    ro.globalenv['pvals'] = ro.FloatVector(list(pvals))
    pval_df['qval'] = ro.r('qvalue(pvals)$qvalues')

    # write to output
    gcs_output = output_path(
        f'fdr_qvals/{cell_type}_qval.tsv',
        'analysis',
    )
    # arrange by ascending q-value
    pval_df = pval_df.sort_values(by='qval', ascending=True)
    pval_df.to_csv(gcs_output, sep='\t', index=False, header=True)


@click.option(
    '--input-dir', help='GCS path directoy to the input gene-level p-value files'
)
@click.option('--cell-types', help='cell type')
@click.option('--chromosomes', help='chromosomes')
@click.command()
def main(input_dir, cell_types, chromosomes):
    """
    Compute Storey's q-values for gene-level p-values
    """
    for cell_type in cell_types.split(','):
        j = get_batch('compute_storey').new_python_job(
            name=f'compute_storey_{cell_type}'
        )
        j.cpu(0.5).memory('lowmem')
        j.call(compute_storey, input_dir, cell_type, chromosomes)
    get_batch('compute_storey').run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
