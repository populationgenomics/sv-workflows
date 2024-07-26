#!/usr/bin/env python3

"""
This script runs R's meta package to generate pooled effect sizes for each eQTL.
Assumes associaTR was run previously on both cohorts and gene lists were generated for each cell type and chromosome.
Outputs a TSV file with the meta-analysis results for each gene.

analysis-runner --dataset "bioheart" --description "meta results runner" --access-level "test" \
    --output-dir "str/associatr/gwas-cell-prop/tob_n1055_and_bioheart_n990" \
     --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    meta_runner.py --results-dir-1=gs://cpg-bioheart-test-analysis/str/associatr/gwas-cell-prop/tob_n1055/results \
    --results-dir-2=gs://cpg-bioheart-test-analysis/str/associatr/gwas-cell-prop/bioheart_n990/results
"""

import click

from cpg_utils.hail_batch import get_batch, image_path, output_path


def run_meta_gen(input_dir_1, input_dir_2, pheno):
    """
    Run meta-analysis (metagen) using R's meta package
    """

    # load packages
    import pandas as pd

    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri


    ro.r('library(meta)')
    ro.r('library(tidyverse)')

    # read in raw associaTR results for each cohort for a particular pheno
    d1 = pd.read_csv(f'{input_dir_1}/{pheno}_fmed_strs.tsv', sep='\t')
    d2 = pd.read_csv(f'{input_dir_2}/{pheno}_fmed_strs.tsv', sep='\t')

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
        f'r_df1 = r_df1 %>% select(chrom, pos, coeff_{pheno}_proportion,se_{pheno}_proportion, motif, period, ref_len,n_samples_tested, `regression_R^2`, allele_frequency)',
    )
    ro.r(
        f'r_df2 = r_df2 %>% select(chrom, pos, coeff_{pheno}_proportion,se_{pheno}_proportion,motif, period, ref_len,n_samples_tested, `regression_R^2`, allele_frequency)',
    )
    ro.r(f'r_df1 = r_df1 %>% rename(coeff_1 =coeff_{pheno}_proportion, se_1 =se_{pheno}_proportion  )')
    ro.r(f'r_df2 = r_df2 %>% rename(coeff_2 =coeff_{pheno}_proportion, se_2 =se_{pheno}_proportion )')
    ro.r('df <- merge(r_df1, r_df2, by = c("chrom", "pos", "motif", "period", "ref_len"))')

    # initialise empty df to store meta results
    ro.r(
        '''
    meta_df <- data.frame(
    chr = character(),
    pos = numeric(),
    n_samples_tested_1 = numeric(),
    n_samples_tested_2 = numeric(),
    coeff_meta = numeric(),
    se_meta = numeric(),
    pval_q_meta = numeric(),
    pval_meta = numeric(),
    lowerCI_meta = numeric(),
    upperCI_meta = numeric(),
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

    new_entry <- data.frame(
        chr = df[i, "chrom"],
        pos = df[i, "pos"],
        n_samples_tested_1 = df[i, "n_samples_tested.x"],
        n_samples_tested_2 = df[i, "n_samples_tested.y"],
        coeff_meta = m.gen$TE.random,
        se_meta = m.gen$seTE.random,
        pval_q_meta = m.gen$pval.Q,
        pval_meta = m.gen$pval.random,
        lowerCI_meta = m.gen$lower.random,
        upperCI_meta = m.gen$upper.random,
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
        f'{output_path(f"meta_results/{pheno}_fmed_strs_meta_results.tsv", "analysis")}',
        sep='\t',
        index=False,
    )


@click.option('--results-dir-1', help='GCS path directory to the raw associatr results for cohort 1')
@click.option('--results-dir-2', help='GCS path directory to the raw associatr results for cohort 2')
@click.command()
def main(
    results_dir_1,
    results_dir_2,
):
    """
    Compute meta-analysis pooled results for each pheno
    """

    for macrogroup in [
        'B',
        'CD4T',
        'CD8T',
        'Cycling',
        'DC',
        'MAIT',
        'Mono',
        'NK',
        'Platelet',
        'Treg',
        'gdT',
        'nan',
        'Eryth',
    ]:
        j = get_batch(name='compute_meta').new_python_job(name=f'compute_meta_{macrogroup}')
        j.cpu(0.25)

        j.image(image_path('r-meta'))
        j.call(run_meta_gen, results_dir_1, results_dir_2, macrogroup)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
