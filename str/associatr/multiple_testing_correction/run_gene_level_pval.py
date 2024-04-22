#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter,too-many-arguments,too-many-locals,too-many-nested-blocks
"""

This script computes gene-level p-values using ACAT (the first step of multiple testing correction).
Assumed input files follow the format of output TSV files from associaTR.
Output is multiple gene-specific TSV files with the gene name in the first column, and gene-level p-value in the second column.
The attributes of the locus with the lowest raw p-value are also stored in the TSV file (coordinates, beta, se, raw pval,r2, motif, ref_len).

analysis-runner --dataset "bioheart" --description "compute gene level pvals" --access-level "test" \
    --output-dir "str/associatr/tob_n1055/results" \
    run_gene_level_pval.py --input-dir=gs://cpg-bioheart-test/str/associatr/tob_n1055/results/v1 \
    --cell-types=CD4_TCM --chromosomes=1 --acat
"""
import click
import numpy as np
import pandas as pd

import hailtop.batch as hb

from cpg_utils.hail_batch import get_batch, image_path

# store a mapping of the key description to the index
VALUES_TO_INDEXES = [
    ('chr', 0),
    ('pos', 1),
    ('n_samples_tested', 3),
    ('raw_pval', 5),
    ('coeff', 6),
    ('se', 7),
    ('r2', 8),
    ('motif', 9),
    ('ref_len', 11),
    ('allele_frequency', 12),
]


def cct(
    gene_name,
    pvals,
    cell_type,
    chromosome,
    row_dict,
    weights=None,
):
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
    from scipy.stats import cauchy
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path


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

    # write to output
    gcs_output = output_path(
        f'gene_level_pvals/acat/{cell_type}/chr{chromosome}/{gene_name}_gene_level_pval.tsv',
        'analysis',
    )
    with to_path(gcs_output).open('w') as f:
        f.write(
            'gene_name\tgene_level_pval\tchr\tpos\tn_samples_tested\tlowest_raw_pval\tcoeff\tse\tr2\tmotif\tref_len\tallele_freq\n',
        )
        f.write(f'{gene_name}\t{pval}\t')
        f.write('\t'.join([str(row_dict[key]) for key, _value in VALUES_TO_INDEXES]) + '\n')


def bonferroni_compute(
    gene_name,
    pvals,
    cell_type,
    chromosome,
    row_dict,
):
    """
    Computes Bonferroni adjusted p-value of the lowest raw p-value for a gene
    """
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    pval = min(pvals) * len(pvals)
    # write to output
    gcs_output = output_path(
        f'gene_level_pvals/bonferroni/{cell_type}/chr{chromosome}/{gene_name}_gene_level_pval.tsv',
        'analysis',
    )
    with to_path(gcs_output).open('w') as f:
        f.write(
            'gene_name\tgene_level_pval\tchr\tpos\tn_samples_tested\tlowest_raw_pval\tcoeff\tse\tr2\tmotif\tref_len\tallele_freq\n',
        )
        f.write(f'{gene_name}\t{pval}\t')
        f.write('\t'.join([str(row_dict[key]) for key, _value in VALUES_TO_INDEXES]) + '\n')


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
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=500,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
@click.option(
    '--acat',
    is_flag=True,
    help='Run ACAT method',
)
@click.option(
    '--bonferroni',
    is_flag=True,
    help='Run Bonferroni method',
)
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
        for chromosome in chromosomes.split(','):
            b = get_batch()
            b.image(image_path('scanpy'))
            gene_files = list(to_path(f'{input_dir}/{cell_type}/chr{chromosome}').glob('*.tsv'))
            for i in range(0, len(gene_files), genes_per_job):
                batch_gene_files = gene_files[i : i + genes_per_job]
                j = b.new_python_job(
                    name=f'Compute gene-level p-values for genes {i+1}-{i+genes_per_job}',
                )
                j.image(image_path('scanpy'))

                j.cpu(0.25).memory('lowmem')
                f = b.new_python_job(
                    name=f'Compute gene-level Bonferroni p-values for genes {i+1}-{i+genes_per_job}',
                )
                f.image(image_path('scanpy'))

                f.cpu(0.25).memory('lowmem')
                for gene_file in batch_gene_files:
                    # read the raw results
                    gene_results = pd.read_csv(gene_file, sep='\t')
                    pvals = gene_results.iloc[:, 5]  # stored in the 6th column
                    # Find and store the attributes of the locus with lowest raw pval
                    # Find the minimum value in column 6
                    min_value = gene_results.iloc[:, 5].min()
                    # Find the rows with the minimum value in column 6
                    min_rows = gene_results[gene_results.iloc[:, 5] == min_value]

                    # create a dictionary of {key: list}
                    row_dict = {key: [] for key, value in VALUES_TO_INDEXES}

                    # populate the dict
                    for _index, row in min_rows.iterrows():
                        for key, value in VALUES_TO_INDEXES:
                            row_dict[key].append(row.iloc[value])

                    pvals = np.array(pvals)
                    gene_name = gene_results.columns[5].split('_')[-1]
                    if acat:
                        j.call(
                            cct,
                            gene_name,
                            pvals,
                            cell_type,
                            chromosome,
                            row_dict,
                        )
                    if bonferroni:
                        f.call(
                            bonferroni_compute,
                            gene_name,
                            pvals,
                            cell_type,
                            chromosome,
                            row_dict,
                        )
                manage_concurrency_for_job(j)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments,too-many-locals,too-many-nested-blocks