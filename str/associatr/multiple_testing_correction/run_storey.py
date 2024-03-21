#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""

This script performs FDR (across gene) correction (the second and final step of multiple testing correction).
Ensure that `run_gene_level_pval.py` has been run to generate the gene-level p-values first.
Output is one TSV file per cell type with three columns - gene name, gene-level p-value (ACAT/Bonferroni), and q-value.

analysis-runner --dataset "bioheart" --description "compute qvals" --access-level "test" \
    --output-dir "str/associatr/rna_pc_calibration/5_pcs/results" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images-dev/r-qvalue:1.0 \
    run_storey.py --input-dir=gs://cpg-bioheart-test/str/associatr/rna_pc_calibration/5_pcs/results/gene_level_pvals \
    --cell-types=CD8_TEM --chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
    --gene-level-correction=bonferroni

"""

import pandas as pd

import click
from cpg_utils.hail_batch import get_batch, output_path
from cpg_utils import to_path
import rpy2.robjects as ro


def compute_storey(input_dir, cell_type, chromosomes, gene_level_correction):
    """
    Compute Storey's q-values for gene-level p-values
    """
    ro.r('library(qvalue)')

    first_iteration = True
    for chromosome in chromosomes.split(','):
        # read in gene-level p-values
        gene_pval_files = list(
            to_path(
                f'{input_dir}/{gene_level_correction}/{cell_type}/chr{chromosome}'
            ).glob('*.tsv')
        )
        if first_iteration:
            pval_df = pd.read_csv(gene_pval_files[0], sep='\t')
            first_iteration = False
        for gene_pval_file in gene_pval_files[1:]:
            pval_df = pd.concat([pval_df, pd.read_csv(gene_pval_file, sep='\t')])

    pvals = pval_df['gene_level_pval']
    if gene_level_correction == 'bonferroni':
        pvals = list(pvals)
        # set p-values > 1 to 1
        pvals = [1 if num > 1 else num for num in pvals]
    ro.globalenv['pvals'] = ro.FloatVector(list(pvals))
    pval_df['qval'] = ro.r('qvalue(pvals)$qvalues')
    pi_o = ro.r('qvalue(pvals)$pi0')
    print(f'pi0: {pi_o}')

    # write to output
    gcs_output = output_path(
        f'fdr_qvals/using_{gene_level_correction}/{cell_type}_qval.tsv',
        'analysis',
    )
    # arrange by ascending q-value
    pval_df = pval_df.sort_values(by='qval', ascending=True)
    pval_df.to_csv(gcs_output, sep='\t', index=False, header=True)


@click.option(
    '--input-dir', help='GCS path directory to the input gene-level p-value files'
)
@click.option(
    '--gene-level-correction',
    type=click.Choice(['acat', 'bonferroni']),
    help='Choose between "acat" or "bonferroni"',
)
@click.option('--cell-types', help='cell type')
@click.option('--chromosomes', help='chromosomes')
@click.command()
def main(input_dir, cell_types, chromosomes, gene_level_correction):
    """
    Compute Storey's q-values for gene-level p-values
    """
    for cell_type in cell_types.split(','):
        j = get_batch(f'compute_storey {gene_level_correction}').new_python_job(
            name=f'compute_storey_{cell_type}'
        )
        j.cpu(1)
        j.call(compute_storey, input_dir, cell_type, chromosomes, gene_level_correction)
    get_batch('compute_storey').run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
