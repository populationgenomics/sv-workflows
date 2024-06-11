#!/usr/bin/env python3

"""

This script will run SusieR, fine-mapping tool.

Required inputs:
- output from`corr_matrix_maker.py` (ie LD matrix)
- associaTR raw outputs (eSNPs and eSTRs combined), preferably also run with `remove_STR_indels.py`.
analysis-runner --dataset "bioheart" \
    --description "Run susieR for eGenes identified by STR analysis" \
    --access-level "test" \
    --image "australia-southeast1-docker.pkg.dev/cpg-common/images-dev/r-meta:susie" \
    --output-dir "str/associatr/test-only" \
    susie_runner.py \
    --celltypes "ASDC" \
    --chromosomes "chr1" \
    --ld-dir "gs://cpg-bioheart-test/str/associatr/fine_mapping/prep_files/test-str-only" \
    --associatr-dir "gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results" \
    --max-parallel-jobs 100
"""

import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def susie_runner(ld_path, associatr_path, celltype, chrom):

    import rpy2.robjects as ro
    import pandas as pd
    from rpy2.robjects import pandas2ri
    from cpg_utils.hail_batch import output_path


    ro.r('library(coloc)')
    ro.r('library(tidyverse)')

    # load in LD matrix and associatr file into R environment
    ld_matrix = pd.read_csv(ld_path, sep='\t')
    gene = ld_path.split('/')[-1].split('_')[0]
    with (ro.default_converter + pandas2ri.converter).context():
        ld_r = ro.conversion.get_conversion().py2rpy(ld_matrix)
    print('loaded in ld_r')
    associatr_sum_stats = pd.read_csv(associatr_path, sep='\t')
    with (ro.default_converter + pandas2ri.converter).context():
        associatr_r = ro.conversion.get_conversion().py2rpy(associatr_sum_stats)
    print('loaded in associatr_r')
    ro.globalenv['associatr_r'] = associatr_r
    ro.globalenv['ld_r'] = ld_r
    ro.globalenv['gene'] = gene

    ro.r(
        '''
    df_subset <- subset(ld_r, select = -X)
    corr_x = as.matrix(df_subset)
    rownames(corr_x) <- NULL
    colnames(corr_x) <- NULL

    # Create a new variable directly with the desired format
    associatr_r$varid <- paste(associatr_r$chr, associatr_r$pos, sep = ":")
    associatr_r$varid <- paste(associatr_r$varid, associatr_r$motif, sep = "_")

    # Create an index based on the order vector
    index <- match(ld_r$X, associatr_r$varid)

    # Reorder the associatr data frame using the index
    df_ordered <- associatr_r[index, ]

    # Run SusieR
    fitted_rss1 <- susie_rss(bhat = df_ordered$coeff_meta, shat = df_ordered$se_meta, n = df_ordered$n_samples_tested_1[1]+df_ordered$n_samples_tested_2[1], R = corr_x, var_y = 1, L = 10,
                         estimate_residual_variance = TRUE)

    # Append SusieR results to dataframe
    df_ordered$susie_pip = susie_get_pip(fitted_rss1, prune_by_cs = TRUE)

    '''
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        susie_associatr_df = ro.conversion.get_conversion().rpy2py(ro.r('df_ordered'))
    print('converted back to pandas df')

    # write to GCS
    susie_associatr_df.to_csv(
        output_path(f"susie/{celltype}/{chrom}/{gene}_100kb.tsv", 'analysis'),
        sep='\t',
        index=False,
    )


@click.option('--celltypes', help='Cell types comma separated')
@click.option('--chromosomes', help='Chromosomes comma separated')
@click.option('--ld-dir', help='Directory to LD correlation matrices')
@click.option('--associatr-dir', help='Directory to associatr outputs')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs', default=500)
@click.command()
def main(celltypes, chromosomes, ld_dir, associatr_dir, max_parallel_jobs):
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
                associatr_path = f'{associatr_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv'

                susie_job = b.new_python_job(
                    f'SusieR for {chrom}:{gene}: {celltype}',
                )
                susie_job.cpu(0.25)
                susie_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/r-meta:susie')
                susie_job.call(susie_runner, ld_file, associatr_path, celltype, chrom)
                manage_concurrency_for_job(susie_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()
