#!/usr/bin/env python3

"""
This script runs R's meta package to generate pooled effect sizes for each eQTL.
Assumes associaTR was run previously on both cohorts and gene lists were generated for each cell type and chromosome.
Outputs a TSV file with the meta-analysis results for each gene.

analysis-runner --dataset "tenk10k" --description "meta results runner" --access-level "test" \
    --output-dir "str/associatr/final_freeze/meta_fixed/cond_analysis_on_snv/bioheart_n975_tob_n950" \
    meta_runner.py --results-dir-1=gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/meta_fixed/cond_analysis_on_snv/bioheart_n975/results/v1-meta-fixed \
    --results-dir-2=gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/meta_fixed/cond_analysis_on_snv/tob_n950/results/v1-meta-fixed
"""

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, image_path, output_path


def run_meta_gen(input_dir_1, input_dir_2, cell_type, chr, gene):
    """
    Run meta-analysis (metagen) using R's meta package
    """

    # load packages
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(meta)')
    ro.r('library(tidyverse)')

    # read in raw associaTR results for each cohort for a particular gene
    d1 = pd.read_csv(f'{input_dir_1}/{cell_type}/{chr}/{gene}_100000bp.tsv', sep='\t',dtype={'locus_filtered': str})
    d2 = pd.read_csv(f'{input_dir_2}/{cell_type}/{chr}/{gene}_100000bp.tsv', sep='\t',dtype={'locus_filtered': str})

    # remove loci that failed to be tested in either dataset
    d1 = d1[d1['locus_filtered'] == 'False']
    d2 = d2[d2['locus_filtered'] == 'False']

    # convert to R dataframe
    with (ro.default_converter + pandas2ri.converter).context():
        r_from_pd_df1 = ro.conversion.get_conversion().py2rpy(d1)

    with (ro.default_converter + pandas2ri.converter).context():
        r_from_pd_df2 = ro.conversion.get_conversion().py2rpy(d2)

    # initialise R variable in the R environment
    ro.globalenv['r_df1'] = r_from_pd_df1
    ro.globalenv['r_df2'] = r_from_pd_df2

    # clean up dfs and merge
    ro.r(
        f'r_df1 = r_df1 %>% select(chrom, pos, coeff_{cell_type}_{chr}_{gene},se_{cell_type}_{chr}_{gene}, motif, period, ref_len,n_samples_tested, `regression_R^2`, allele_frequency)',
    )
    ro.r(
        f'r_df2 = r_df2 %>% select(chrom, pos, coeff_{cell_type}_{chr}_{gene},se_{cell_type}_{chr}_{gene},motif, period, ref_len,n_samples_tested, `regression_R^2`, allele_frequency)',
    )
    ro.r(f'r_df1 = r_df1 %>% rename(coeff_1 =coeff_{cell_type}_{chr}_{gene}, se_1 =se_{cell_type}_{chr}_{gene}  )')
    ro.r(f'r_df2 = r_df2 %>% rename(coeff_2 =coeff_{cell_type}_{chr}_{gene}, se_2 =se_{cell_type}_{chr}_{gene}  )')
    ro.r('df <- merge(r_df1, r_df2, by = c("chrom", "pos", "motif", "period", "ref_len"))')

    # initialise empty df to store meta results
    ro.r(
        '''
    meta_df <- data.frame(
    chr = character(),
    pos = numeric(),
    n_samples_tested_1 = numeric(),
    n_samples_tested_2 = numeric(),
    coeff_meta_random = numeric(),
    se_meta_random = numeric(),
    pval_q_meta = numeric(),
    q_meta = numeric(),
    tau_meta = numeric(),
    i2_meta = numeric(),
    h_meta = numeric(),
    pval_meta_random = numeric(),
    pval_meta_fixed = numeric(),
    coeff_meta_fixed = numeric(),
    se_meta_fixed = numeric(),
    lowerCI_meta_random = numeric(),
    upperCI_meta_random = numeric(),
    r2_1 = numeric(),
    r2_2 = numeric(),
    motif = character(),
    period = numeric(),
    ref_len = numeric(),
    allele_frequency_1 = numeric(),
    allele_frequency_2 = numeric()
    )
    ''',
    )

    # run meta

    ro.r(
        '''
    for (i in 1:nrow(df)) {
    # Extract values for cohort 1
    row_cohort1 <- data.frame(
        cohort = 1,
        se = df[i, "se_1"],
        coeff = df[i, "coeff_1"]
    )
    # Extract values for cohort 2
    row_cohort2 <- data.frame(
        cohort = 2,
        se = df[i, "se_2"],
        coeff = df[i, "coeff_2"]
    )
    result_df <- rbind(row_cohort1, row_cohort2)

    skip_to_next <- FALSE

    tryCatch({
        m.gen <- metagen(result_df$coeff, result_df$se, random = FALSE, control = list(stepadj = 0.5))
    }, error = function(e) {
        cat("Error at iteration:", i, "\n")
        skip_to_next <<- TRUE
    })

    if (skip_to_next) next

    new_entry = data.frame(
        chr = df[i, "chrom"],
        pos = df[i, "pos"],
        n_samples_tested_1 = df[i, "n_samples_tested.x"],
        n_samples_tested_2 = df[i, "n_samples_tested.y"],
        coeff_meta_random = m.gen$TE.random,
        se_meta_random = m.gen$seTE.random,
        pval_q_meta = m.gen$pval.Q,
        q_meta = m.gen$Q,
        tau_meta = m.gen$tau,
        i2_meta = m.gen$I2,
        h_meta = m.gen$H,
        pval_meta_random = m.gen$pval.random,
        pval_meta_fixed = m.gen$pval.fixed,
        coeff_meta_fixed = m.gen$TE.fixed,
        se_meta_fixed = m.gen$seTE.fixed,
        lowerCI_meta_random = m.gen$lower.random,
        upperCI_meta_random = m.gen$upper.random,
        r2_1 = df[i, "regression_R^2.x"],
        r2_2 = df[i, "regression_R^2.y"],
        motif = df[i, "motif"],
        period = df[i, "period"],
        ref_len = df[i, "ref_len"],
        allele_frequency_1 = df[i, "allele_frequency.x"],
        allele_frequency_2 = df[i, "allele_frequency.y"]
    )

    meta_df <- rbind(meta_df, new_entry)
    }
    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_meta_df = ro.conversion.get_conversion().rpy2py(ro.r('meta_df'))

    # write to GCS
    pd_meta_df.to_csv(
        f'{output_path(f"meta_results/{cell_type}/{chr}/{gene}_100000bp_meta_results.tsv", "analysis")}',
        sep='\t',
        index=False,
    )


@click.option('--results-dir-1', help='GCS path directory to the raw associatr results for cohort 1')
@click.option('--results-dir-2', help='GCS path directory to the raw associatr results for cohort 2')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs to run', default=500)
@click.option('--always-run', is_flag=True, help='Set job to always run')
@click.command()
def main(
    results_dir_1,
    results_dir_2,
    max_parallel_jobs,
    always_run,
):
    """
    Compute meta-analysis pooled results for each gene
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

    celltypes = 'gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive'
    celltypes = celltypes.split(',')
    for cell_type in celltypes:
        gene_files = list(map(str, to_path(f'{results_dir_1}/{cell_type}').rglob('*.tsv')))
        for gene_file in gene_files:
            gene = gene_file.split('/')[-1].split('_')[0]
            chrom = gene_file.split('/')[-2]
            if to_path(
                output_path(
                    f"meta_results/{cell_type}/{chrom}/{gene}_100000bp_meta_results.tsv",
                    "analysis",
                ),
            ).exists():
                continue
            j = get_batch(name='compute_meta').new_python_job(name=f'compute_meta_{cell_type}_{chrom}_{gene}')
            j.cpu(0.25)
            if always_run:
                j.always_run()
            j.image(image_path('r-meta'))
            j.call(run_meta_gen, results_dir_1, results_dir_2, cell_type, chrom, gene)
            manage_concurrency_for_job(j)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
