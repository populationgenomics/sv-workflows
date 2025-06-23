#!/usr/bin/env python3

"""

This script will run SusieR, fine-mapping tool.

Required inputs:
- output from`corr_matrix_maker.py` (ie LD matrix)
- associaTR raw outputs (eSNPs and eSTRs combined), preferably also run with `remove_STR_indels.py`.
analysis-runner --dataset "tenk10k" \
    --description "Run susieR for eGenes identified by STR analysis" \
    --access-level "test" \
    --image "australia-southeast1-docker.pkg.dev/cpg-common/images/r-meta:susie" \
    --output-dir "str/associatr/final_freeze/fine_mapping" \
    susie_runner.py \
    --celltypes "CD4_TCM" \
    --chromosomes "chr1,chr2,chr3,chr4,chr5" \
    --ld-dir "gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/prep_files/v2/correlation_matrix" \
    --associatr-dir "gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/trs_snps/rm_str_indels_dup_strs_v3" \
    --max-parallel-jobs 100 --num-causal-variants=1
"""

import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path


def susie_runner(ld_path, associatr_path, celltype, chrom, num_iterations, num_causal_variants):
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(susieR)')
    ro.r('library(tidyverse)')

    # load in LD matrix and associatr file into R environment
    # Load LD matrix
    ld_matrix = pd.read_csv(ld_path, sep='\t')
    ld_matrix = ld_matrix.rename(columns={'Unnamed: 0': 'X'})
    gene = ld_path.split('/')[-1].split('_')[0]
    # Load and clean sumstats
    associatr_sum_stats = pd.read_csv(associatr_path, sep='\t')
    associatr_sum_stats = associatr_sum_stats.dropna(subset=['coeff_meta_fixed', 'se_meta_fixed'])

    # Generate varid for matching
    associatr_sum_stats['varid'] = associatr_sum_stats['chr'].astype(str) + ":" + associatr_sum_stats['pos'].astype(str)
    associatr_sum_stats['varid'] = associatr_sum_stats['varid'] + "_" + associatr_sum_stats['motif']

    # Filter LD matrix to match only valid rows
    ld_matrix = ld_matrix[ld_matrix['X'].isin(associatr_sum_stats['varid'])]

    # Subset associatr to only what matches LD
    associatr_sum_stats = associatr_sum_stats[associatr_sum_stats['varid'].isin(ld_matrix['X'])]

    # Reorder associatr to match LD
    associatr_sum_stats = associatr_sum_stats.set_index('varid').loc[ld_matrix['X']].reset_index()

    # Drop 'X' column from LD matrix to make it numeric
    corr_x = ld_matrix.drop(columns='X')

    # Convert to R
    with (ro.default_converter + pandas2ri.converter).context():
        ro.globalenv['ld_r'] = ro.conversion.get_conversion().py2rpy(corr_x)
        ro.globalenv['associatr_r'] = ro.conversion.get_conversion().py2rpy(associatr_sum_stats)
        ro.globalenv['gene'] = gene
        ro.globalenv['num_iterations'] = num_iterations
        ro.globalenv['num_causal_variants'] = num_causal_variants

    # R code stays clean now
    ro.r('''
    corr_x <- as.matrix(ld_r)
    rownames(corr_x) <- NULL
    colnames(corr_x) <- NULL

    # Run SuSiE
    fitted_rss1 <- susie_rss(
        bhat = associatr_r$coeff_meta_fixed,
        shat = associatr_r$se_meta_fixed,
        n = associatr_r$n_samples_tested_1[1] + associatr_r$n_samples_tested_2[1],
        R = corr_x,
        var_y = 1,
        L = num_causal_variants,
        max_iter = num_iterations
    )
   # Append SusieR results to dataframe
    df_ordered$susie_pip = susie_get_pip(fitted_rss1, prune_by_cs = TRUE)

    # Capture the raw outputs
    raw_output = capture.output(fitted_rss1)

    # Capture the summary outputs (try catch, in case it is empty)
    p4 <- tryCatch({

    summary_df <- as.data.frame(as.matrix(summary(fitted_rss1)))
    summary_df <- summary_df[order(summary_df$variable), ]
    summary_df <- summary_df[, -1]
    df_ordered <- cbind(df_ordered, summary_df)
    }, error = function(e) {

    print(paste("An error occurred:", e))
    0
    })


    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        susie_associatr_df = ro.conversion.get_conversion().rpy2py(ro.r('df_ordered'))
    print('converted back to pandas df')

    # convert raw output to python
    raw_output_python = ro.r('raw_output')

    # write raw output to GCS
    with to_path(output_path(f"susie/{celltype}/{chrom}/{gene}_100kb_output.txt", 'analysis')).open('w') as file:
        file.write(str(raw_output_python))

    # write dataframe to GCS
    susie_associatr_df.to_csv(
        output_path(f"susie_summstats/{celltype}/{chrom}/{gene}_100kb.tsv", 'analysis'),
        sep='\t',
        index=False,
    )


@click.option('--celltypes', help='Cell types comma separated')
@click.option('--chromosomes', help='Chromosomes comma separated')
@click.option('--ld-dir', help='Directory to LD correlation matrices')
@click.option('--associatr-dir', help='Directory to associatr outputs')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs', default=500)
@click.option('--num_iterations', help='Number of iterations for SusieR', default=100)
@click.option('--susie-cpu', help='CPU for SusieR job', default=0.25)
@click.option('--num-causal-variants', help='Number of causal variants to estimate', default=10)
@click.option('--always-run', help='Job set to always run', is_flag=True)
@click.command()
def main(
    celltypes,
    chromosomes,
    ld_dir,
    associatr_dir,
    max_parallel_jobs,
    num_iterations,
    susie_cpu,
    num_causal_variants,
    always_run,
):
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    b = get_batch(name='Run susieR')

    for celltype in celltypes.split(','):
        for chrom in chromosomes.split(','):
            ld_files = list(to_path(f'{ld_dir}/{celltype}/{chrom}').glob('*.tsv'))
            for ld_file in ld_files:  # for each gene (each has its own LD file)
                ld_file = str(ld_file)
                gene = ld_file.split('/')[-1].split('_')[0]
                print(f'Processing {gene}...')
                if to_path(
                    output_path(f"susie_summstats/{celltype}/{chrom}/{gene}_100kb.tsv", 'analysis'),
                ).exists():
                    continue
                associatr_path = f'{associatr_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv'

                susie_job = b.new_python_job(
                    f'SusieR for {chrom}:{gene}: {celltype}',
                )
                susie_job.cpu(susie_cpu)
                if always_run:
                    susie_job.always_run()
                susie_job.call(
                    susie_runner,
                    ld_file,
                    associatr_path,
                    celltype,
                    chrom,
                    num_iterations,
                    num_causal_variants,
                )
                manage_concurrency_for_job(susie_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()
