#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""

This script computes gene-level p-values using ACAT (the first step of multiple testing correction).
Assumed input files follow the format of output TSV files from meta_runner.py (ie post meta-analysis output files).
Output is multiple gene-specific TSV files with the gene name in the first column, and gene-level p-value in the second column.
The attributes of the locus with the lowest raw p-value are also stored in the TSV file (coordinates, pooled _beta, pooled_se, pooled_pval,pooled_pval_q, motif, ref_len).

analysis-runner --dataset "tenk10k" --description "compute gene level pvals" --access-level "test" \
    --output-dir "str/associatr/rna_calib/tob_n950_and_bioheart_n975/pc1" \
    run_gene_level_pvals_meta.py --input-dir=gs://cpg-tenk10k-test-analysis/str/associatr/rna_calib/tob_n950_and_bioheart_n975/pc1/meta_results \
    --cell-types=CD14_Mono,CD4_TCM,CD8_TEM \
    --chromosomes=1 --acat
"""
import logging
from copy import deepcopy

import click
import numpy as np

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, reset_batch

# store a mapping of the key description to the index
VALUES_TO_INDEXES = [
    ('chr', 0),
    ('pos', 1),
    ('n_samples_tested_1', 2),
    ('n_samples_tested_2', 3),
    ('coeff_meta_fixed', 13),
    ('se_meta_fixed', 14),
    ('pval_q', 6),
    ('pval_meta_fixed',12),
    ('r2_1', 17),
    ('r2_2', 18),
    ('motif', 19),
    ('ref_len', 21),
    ('allele_frequency_1', 22),
    ('allele_frequency_2', 23),
]


def process_single_file(gene_file: str):
    """
    reads a file, and pulls out all the good stuff

    Args:
        gene_file (str):

    Returns:
        pvals, gene_name, row_dict
    """
    import pandas as pd

    # read the raw results
    gene_results = pd.read_csv(gene_file, sep='\t')
    pvals = gene_results.iloc[:, 12]  # stored in the 13th column
    # Find and store the attributes of the locus with lowest raw pval
    # Find the minimum value in column 13
    min_value = gene_results.iloc[:, 12].min()
    # Find the rows with the minimum value in column 13
    min_rows = gene_results[gene_results.iloc[:, 12] == min_value]

    # create a dictionary of {key: list}
    row_dict: dict[str, list] = {key: [] for key, value in VALUES_TO_INDEXES}

    # populate the dict
    for _index, row in min_rows.iterrows():
        for key, value in VALUES_TO_INDEXES:
            row_dict[key].append(row.iloc[value])

    pvals = np.array(pvals)
    gene_name = gene_file.split('/')[-1].split('_')[0]
    return pvals, gene_name, row_dict


def cct(gene_files: list[str], cell_type: str, chromosome: str, og_weights=None):
    """
    takes a list of gene files to process in this job
    moves processing out of the driver job

    weights were breaking - this script never supplies a value to weights, but
    the looping caused weights to be taken from the previous round's value

    Args:
        gene_files (list): the list of files to process
        cell_type (str):
        chromosome (str):
        og_weights (): None

    Returns:

    """
    """
    Code adapted from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value
    #' An analytical p-value combination method using the Cauchy distribution
    #'
    #' The code{CCT} function takes in a numeric vector of p-values, a numeric
    #' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
    #' @param pvals a numeric vector of p-values, where each of the element is
    #' between 0 to 1, to be combined.
    #' @param weights a numeric vector of non-negative weights. If code{NULL}, the
    #' equal weights are assumed.
    #' @return the aggregated p-value combining p-values from the vector code{pvals}.
    #' @examples pvalues <- c(2e-02,4e-04,0.2,0.1,0.8)
    #' @examples CCT(pvals=pvalues)
    #' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
    #' with analytic p-value calculation under arbitrary dependency structures.
    #' emph{Journal of the American Statistical Association 115}(529), 393-402.
    #' (href{https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485}{pub})
    #' @export
    R code is implemented in python
    """
    # Import here as a PythonJob's function must stand alone
    from scipy.stats import cauchy  # noqa: PLC0415

    from cpg_utils import to_path
    from cpg_utils.config import output_path

    for gene_file in gene_files:
        pvals, gene_name, row_dict = process_single_file(gene_file)

        # specify output path
        gcs_output = to_path(
            output_path(
                f'gene_level_pvals/acat/{cell_type}/chr{chromosome}/{gene_name}_gene_level_pval.tsv',
                'analysis',
            ),
        )
        if gcs_output.exists():
            print(f'{gene_file} already processed. Skipping...')
            continue
        weights = deepcopy(og_weights)
        print(f'processing {gene_file}')

        # remove NA values - associaTR reports pval as NA if locus was thrown out (not tested)
        pvals = pvals[~np.isnan(pvals)]

        # check if all p-values are between 0 and 1
        if ((pvals < 0).sum() + (pvals > 1).sum()) > 0:
            raise ValueError('All p-values must be between 0 and 1!')

        # check if there are p-values that are either exactly 0 or 1.
        is_zero = (pvals == 0).sum() >= 1
        is_one = (pvals == 1).sum() >= 1

        # check the validity of weights (default: equal weights) and standardize them.
        if weights is None:
            weights = np.repeat(1 / len(pvals), len(pvals))
        elif len(weights) != len(pvals):
            raise ValueError('The length of weights should be the same as that of the p-values!')
        elif (weights < 0).sum() > 0:
            raise ValueError('All the weights must be positive!')
        else:
            weights = weights / np.sum(weights)

        # check if there are very small non-zero p-values
        is_small = pvals < 1e-16
        if is_small.sum() == 0:
            cct_stat = np.sum(weights * np.tan((0.5 - pvals) * np.pi))
        else:
            cct_stat = np.sum((weights[is_small] / pvals[is_small]) / np.pi)
            cct_stat += np.sum(weights[~is_small] * np.tan((0.5 - pvals[~is_small]) * np.pi))

        if is_zero:
            pval = 0
        elif is_one:
            print('There are p-values exactly equal to 1!')
            pval = min(1, min(pvals) * len(pvals))

        # check if the test statistic is very large.
        elif cct_stat > 1e15:
            pval = (1 / cct_stat) / np.pi
        else:
            pval = 1 - cauchy.cdf(cct_stat)

        with to_path(gcs_output).open('w') as f:
            f.write(
                'gene_name\tgene_level_pval\tchr\tpos\tn_samples_tested_1\tn_samples_tested_2\tcoeff_meta_fixed\tse_meta_fixed\tpval_q\tpval_meta_fixed\tr2_1\tr2_2\tmotif\tref_len\tallele_freq_1\tallele_freq_2\n',
            )

            f.write(f'{gene_name}\t{pval}\t')
            f.write('\t'.join([str(row_dict[key]) for key, _value in VALUES_TO_INDEXES]) + '\n')


def bonferroni_compute(gene_files, cell_type, chromosome):
    """
    Computes Bonferroni adjusted p-value of the lowest raw p-value for a gene
    """
    import pandas as pd

    from cpg_utils import to_path
    from cpg_utils.config import output_path

    for gene_file in gene_files:
        # read the raw results from a file
        pvals, gene_name, row_dict = process_single_file(gene_file)

        # write to output
        gcs_output = to_path(
            output_path(
                f'gene_level_pvals/bonferroni/{cell_type}/chr{chromosome}/{gene_name}_gene_level_pval.tsv',
                'analysis',
            ),
        )
        if gcs_output.exists():
            print(f'{gene_file} already processed. Skipping...')
            continue

        pval = min(pvals) * len(pvals)

        with to_path(gcs_output).open('w') as f:
            f.write(
                'gene_name\tgene_level_pval\tchr\tpos\tn_samples_tested_1\tn_samples_tested_2\tcoeff_meta_fixed\tse_meta_fixed\tpval_q\tpval_meta_fixed\tr2_1\tr2_2\tmotif\tref_len\tallele_freq_1\tallele_freq_2\n',
            )
            f.write(f'{gene_name}\t{pval}\t')
            f.write('\t'.join([str(row_dict[key]) for key, _value in VALUES_TO_INDEXES]) + '\n')


@click.option('--input-dir', help='GCS path to the raw results of associaTR')
@click.option('--cell-types', help='Name of the cell type, comma separated if multiple')
@click.option(
    '--chromosomes',
    help='Chromosome number eg 1, comma separated if multiple',
)
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=18,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
@click.option('--acat', is_flag=True, help='Run ACAT method')
@click.option('--bonferroni', is_flag=True, help='Run Bonferroni method')
@click.command()
def main(input_dir, cell_types, chromosomes, max_parallel_jobs, acat, bonferroni):
    """
    Compute gene-level p-values
    """
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    genes_per_job = 70

    for cell_type in cell_types.split(','):
        # run each cell-type batch completely separately to help with job scaling
        _dependent_jobs = []
        reset_batch()
        b = get_batch()

        for chromosome in chromosomes.split(','):
            gene_files = list(map(str, to_path(f'{input_dir}/{cell_type}/chr{chromosome}').glob('*.tsv')))

            # split the list of files into chunks
            for i in range(0, len(gene_files), genes_per_job):
                batch_gene_files = gene_files[i : i + genes_per_job]

                if acat:
                    j = b.new_python_job(
                        name=f'Compute gene-level p-values for genes {i+1}-{i+genes_per_job} in {cell_type}:{chromosome}',
                    )
                    j.cpu(0.25).memory('lowmem')
                    j.call(cct, batch_gene_files, cell_type, chromosome)
                    manage_concurrency_for_job(j)

                if bonferroni:
                    f = b.new_python_job(
                        name=f'Compute gene-level Bonferroni p-values for genes {i+1}-{i+genes_per_job} in {cell_type}:{chromosome}',
                    )
                    f.cpu(0.25).memory('lowmem')
                    f.call(bonferroni_compute, batch_gene_files, cell_type, chromosome)
                    manage_concurrency_for_job(f)
        b.run(wait=False)


if __name__ == '__main__':
    # catch the logging emitted by batch generation
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
