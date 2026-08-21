#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
Compute gene-level p-values from the SV meta-analysis outputs (first step of multiple testing
correction). SV counterpart of str/associatr/meta_analysis/run_gene_level_pvals_meta.py.

Input files are the per-gene TSVs written by meta_runner_sv.py. Output is one TSV per gene: the
gene name, the gene-level p-value, and the attributes of the locus with the lowest pooled p-value.

Two differences from the STR version, both forced by the SV output format:

1. Columns are located by NAME, not by position. The STR script hardcodes indices (pval at
   iloc[:, 12]); in the SV meta output pval_meta_fixed sits at index 37, and the layout has
   already changed once during development. Name-based access also means adding or removing a
   column upstream cannot silently shift the wrong values into the output.

2. The carried locus attributes are SV-appropriate. There is no motif/period/ref_len to report;
   instead each row identifies the federated variant (federate_id, the two cohort VARIDs, their
   positions and SVTYPEs) plus pair_class, so a gene-level hit can be traced straight back to the
   pairing table.

The ACAT/CCT implementation is copied unchanged from the STR script (a port of STAAR's CCT with
their documented modification: when an individual p-value is exactly 1, fall back to the minimum
p-value), so gene-level p-values stay methodologically comparable across the STR and SV analyses.

analysis-runner --dataset "tenk10k" --description "compute SV gene level pvals" --access-level "test" \
    --output-dir "tenk10k-sv/sv/meta_analysis/sensitivity_analysis/bioheart_n968_and_tob_n935/10pc/v1" \
    python3 run_gene_level_pvals_meta_sv.py \
    --input-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/sensitivity_analysis/bioheart_n968_and_tob_n935/10pc/v1/meta_results \
    --cell-types=CD8_TEM \
    --chromosomes=chr1 --acat
"""

import logging
from copy import deepcopy

import click
import numpy as np

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, reset_batch

# Attributes of the lowest-pooled-p locus carried into the gene-level output, by column NAME.
# The output header is generated from this tuple, so the two cannot drift apart.
CARRY_COLUMNS = (
    'federate_id',
    'pair_class',
    'chr',
    'bioheart_vid',
    'bioheart_pos',
    'bioheart_svtype',
    'tob_vid',
    'tob_pos',
    'tob_svtype',
    'n_samples_tested_1',
    'n_samples_tested_2',
    'coeff_meta_fixed',
    'se_meta_fixed',
    'pval_meta_fixed',
    'pval_q_meta',
    'i2_meta',
    'r2_1',
    'r2_2',
    'allele_frequency_1',
    'allele_frequency_2',
    'abs_alt_af_diff',
)

OUTPUT_HEADER = '\t'.join(('gene_name', 'gene_level_pval', 'n_loci') + CARRY_COLUMNS) + '\n'


def process_single_file(gene_file: str, pval_column: str = 'pval_meta_fixed', exclude_dup_cnv: bool = False):
    """Read one meta-analysis gene TSV and pull out the p-values plus the top locus.

    Args:
        gene_file (str): path to a *_meta_results.tsv written by meta_runner_sv.py
        pval_column (str): which pooled p-value to combine ('pval_meta_fixed' or 'pval_meta_random')
        exclude_dup_cnv (bool): drop DUP-vs-CNV pairs before combining (sensitivity check)

    Returns:
        pvals, gene_name, row_dict
    """
    import pandas as pd  # noqa: PLC0415

    gene_results = pd.read_csv(gene_file, sep='\t')

    missing = [c for c in (pval_column,) + CARRY_COLUMNS if c not in gene_results.columns]
    if missing:
        raise ValueError(f'{gene_file} is missing expected column(s): {missing}')

    if exclude_dup_cnv:
        gene_results = gene_results[gene_results['pair_class'] != 'dup_cnv']

    pvals = np.array(gene_results[pval_column])

    # attributes of the locus with the lowest pooled p-value; ties keep every tied row
    row_dict: dict[str, list] = {key: [] for key in CARRY_COLUMNS}
    if len(gene_results):
        min_rows = gene_results[gene_results[pval_column] == gene_results[pval_column].min()]
        for _index, row in min_rows.iterrows():
            for key in CARRY_COLUMNS:
                row_dict[key].append(row[key])

    gene_name = gene_file.split('/')[-1].split('_')[0]
    return pvals, gene_name, row_dict


def _format_row(gene_name, pval, n_loci, row_dict):
    """Gene-level result as one TSV line. Tied top loci are comma-joined within each field."""
    values = [','.join(str(v) for v in row_dict[key]) for key in CARRY_COLUMNS]
    return f'{gene_name}\t{pval}\t{n_loci}\t' + '\t'.join(values) + '\n'


def cct(  # noqa: C901
    gene_files: list[str],
    cell_type: str,
    chromosome: str,
    pval_column: str = 'pval_meta_fixed',
    exclude_dup_cnv: bool = False,
    og_weights=None,
):
    """
    takes a list of gene files to process in this job
    moves processing out of the driver job

    weights were breaking - this script never supplies a value to weights, but
    the looping caused weights to be taken from the previous round's value
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
    #' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
    #' with analytic p-value calculation under arbitrary dependency structures.
    #' emph{Journal of the American Statistical Association 115}(529), 393-402.
    #' @export
    R code is implemented in python
    """
    # Import here as a PythonJob's function must stand alone
    from scipy.stats import cauchy  # noqa: PLC0415

    from cpg_utils import to_path  # noqa: PLC0415
    from cpg_utils.config import output_path  # noqa: PLC0415

    for gene_file in gene_files:
        pvals, gene_name, row_dict = process_single_file(gene_file, pval_column, exclude_dup_cnv)

        gcs_output = to_path(
            output_path(
                f'gene_level_pvals/acat/{cell_type}/{chromosome}/{gene_name}_gene_level_pval.tsv',
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

        # nothing left to combine (all NA, or every locus dropped by --exclude-dup-cnv)
        if len(pvals) == 0:
            print(f'  no usable p-values in {gene_file}; skipping')
            continue

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
            f.write(OUTPUT_HEADER)
            f.write(_format_row(gene_name, pval, len(pvals), row_dict))


def bonferroni_compute(
    gene_files,
    cell_type,
    chromosome,
    pval_column: str = 'pval_meta_fixed',
    exclude_dup_cnv: bool = False,
):
    """Bonferroni-adjust the lowest pooled p-value for each gene."""
    from cpg_utils import to_path  # noqa: PLC0415
    from cpg_utils.config import output_path  # noqa: PLC0415

    for gene_file in gene_files:
        pvals, gene_name, row_dict = process_single_file(gene_file, pval_column, exclude_dup_cnv)

        gcs_output = to_path(
            output_path(
                f'gene_level_pvals/bonferroni/{cell_type}/{chromosome}/{gene_name}_gene_level_pval.tsv',
                'analysis',
            ),
        )
        if gcs_output.exists():
            print(f'{gene_file} already processed. Skipping...')
            continue

        pvals = pvals[~np.isnan(pvals)]
        if len(pvals) == 0:
            print(f'  no usable p-values in {gene_file}; skipping')
            continue

        pval = min(min(pvals) * len(pvals), 1.0)

        with to_path(gcs_output).open('w') as f:
            f.write(OUTPUT_HEADER)
            f.write(_format_row(gene_name, pval, len(pvals), row_dict))


@click.option('--input-dir', required=True, help='GCS path to the meta_results dir from meta_runner_sv.py')
@click.option('--cell-types', required=True, help='Name of the cell type, comma separated if multiple')
@click.option('--chromosomes', required=True, help='Chromosome eg chr1 or 1, comma separated if multiple')
@click.option(
    '--pval-column',
    default='pval_meta_fixed',
    type=click.Choice(['pval_meta_fixed', 'pval_meta_random']),
    help='which pooled p-value to combine',
)
@click.option(
    '--exclude-dup-cnv',
    is_flag=True,
    help='drop DUP-vs-CNV pairs before combining (sensitivity check on the biallelic/multiallelic coding)',
)
@click.option('--genes-per-job', type=int, default=70, help='genes processed per job')
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=18,
    help='To avoid exceeding Google Cloud quotas, set this concurrency as a limit.',
)
@click.option('--acat', is_flag=True, help='Run ACAT method')
@click.option('--bonferroni', is_flag=True, help='Run Bonferroni method')
@click.command()
def main(
    input_dir,
    cell_types,
    chromosomes,
    pval_column,
    exclude_dup_cnv,
    genes_per_job,
    max_parallel_jobs,
    acat,
    bonferroni,
):
    """Compute gene-level p-values from the SV meta-analysis outputs."""
    if not (acat or bonferroni):
        raise click.UsageError('nothing to do: pass --acat and/or --bonferroni')

    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """To avoid having too many jobs running at once, we have to limit concurrency."""
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    for cell_type in cell_types.split(','):
        # run each cell-type batch completely separately to help with job scaling
        _dependent_jobs = []
        reset_batch()
        b = get_batch()

        for chromosome in chromosomes.split(','):
            # accept either 'chr1' or '1'; meta_runner_sv.py writes 'chr1'-style directories
            chrom = chromosome if chromosome.startswith('chr') else f'chr{chromosome}'
            gene_files = sorted(map(str, to_path(f'{input_dir}/{cell_type}/{chrom}').glob('*.tsv')))
            print(f'{cell_type} {chrom}: {len(gene_files):,} gene files')
            if not gene_files:
                continue

            for i in range(0, len(gene_files), genes_per_job):
                batch_gene_files = gene_files[i : i + genes_per_job]
                span = f'{i + 1}-{min(i + genes_per_job, len(gene_files))}'

                if acat:
                    j = b.new_python_job(name=f'ACAT gene-level p-values {span} in {cell_type}:{chrom}')
                    j.cpu(0.25).memory('lowmem')
                    j.call(cct, batch_gene_files, cell_type, chrom, pval_column, exclude_dup_cnv)
                    manage_concurrency_for_job(j)

                if bonferroni:
                    f = b.new_python_job(name=f'Bonferroni gene-level p-values {span} in {cell_type}:{chrom}')
                    f.cpu(0.25).memory('lowmem')
                    f.call(bonferroni_compute, batch_gene_files, cell_type, chrom, pval_column, exclude_dup_cnv)
                    manage_concurrency_for_job(f)

        b.run(wait=False)


if __name__ == '__main__':
    # catch the logging emitted by batch generation
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
