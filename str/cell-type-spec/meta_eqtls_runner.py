#!/usr/bin/env python3

"""
This script performs meta-analysis of eQTLs across cell types to assess cell type-specificity of eQTLs..
Assumes that file_prep.py has been run.


This script runs R's meta package to generate pooled effect sizes for each eQTL tested between pairs of cell types.

analysis-runner --dataset "bioheart" --description "meta_eqtls_runner" --access-level "test" \
--output-dir "str/associatr/cell-type-spec" meta_eqtls_runner.py --file-input-dir=gs://cpg-bioheart-test/str/associatr/cell-type-spec/prep_files \
--cell-types=ASDC

## --cell-types=CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,CD4_TCM_permuted,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC


"""


def run_meta_gen(file_path, cell_type):
    """
    Run meta-analysis (metagen) using R's meta package
    """

    # load packages
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(meta)')
    ro.r('library(tidyverse)')

    # read in raw associaTR results for each cohort for a particular gene
    d1 = pd.read_csv(file_path)

    # convert to R dataframe
    with (ro.default_converter + pandas2ri.converter).context():
        r_from_pd_df1 = ro.conversion.get_conversion().py2rpy(d1)

    # initialise R variable in the R environment
    ro.globalenv['df'] = r_from_pd_df1

    # initialise empty df to store meta results
    ro.r(
        '''
    meta_df <- data.frame(
    chrom = character(),
    pos = numeric(),
    end = numeric()
    gene_name = character(),
    motif = character(),
    celltype_main = character(),
    coeff_main = numeric(),
    se_main = numeric(),
    pval_main = numeric(),
    cell_type2 = character(),
    coeff_2 = numeric(),
    se_2 = numeric(),
    coeff_meta = numeric(),
    se_meta = numeric(),
    pval_q_meta = numeric(),
    pval_meta = numeric(),
    lowerCI_meta = numeric(),
    upperCI_meta = numeric()
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
        se = df[i, "se_main"],
        coeff = df[i, "coeff_main"]
    )
    # Extract values for cohort 2
    row_cohort2 <- data.frame(
        cohort = 2,
        se = df[i, "se_2"],
        coeff = df[i, "coeff_2"]
    )
    result_df <- rbind(row_cohort1, row_cohort2)
    m.gen = metagen(result_df$coeff, result_df$se, random = TRUE)
    print('Meta-analysis ran')


    new_entry = data.frame(
        chrom = df[i, "chrom"],
        pos = df[i, "pos"],
        end = df[i, "end"],
        gene_name = df[i, "gene_name"],
        motif = df[i, "motif"],
        celltype_main = df[i, "celltype_main"],
        coeff_main = df[i, "coeff_main"],
        se_main = df[i, "se_main"],
        pval_main = df[i, "pval_main"],
        cell_type2 = df[i, "cell_type2"],
        coeff_2 = df[i, "coeff_2"],
        se_2 = df[i, "se_2"],
        coeff_meta = m.gen$TE.random,
        se_meta = m.gen$seTE.random,
        pval_q_meta = m.gen$pval.Q,
        pval_meta = m.gen$pval.random,
        lowerCI_meta = m.gen$lower.random,
        upperCI_meta = m.gen$upper.random)

    print('New entry created')

    meta_df = rbind(meta_df, new_entry)
    }
    print('Meta-analysis bound to meta_df')
    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_meta_df = ro.conversion.get_conversion().rpy2py(ro.r('meta_df'))

    # write to GCS
    pd_meta_df.to_csv(
        f'{output_path(f"meta_results/{cell_type}/meta_results.tsv", "analysis")}',
        sep='\t',
        index=False,
    )


import click

from cpg_utils.hail_batch import get_batch,image_path


@click.option('--file-input-dir', required=True, help='Directory containing input files for meta-analysis')
@click.option('--cell-types', required=True, help='Comma-separated list of cell types to compare')
@click.command()
def main(file_input_dir, cell_types):
    b = get_batch(name='meta_eqtl_cell_spec_runner')
    for cell_type in cell_types.split(','):
        file_path = f'{file_input_dir}/{cell_type}/meta_input_df.csv'
        meta_job = b.new_python_job(name=f'{cell_type}_meta_eqtl_cell_spec_runner')
        meta_job.image(image_path('r-meta'))
        meta_job.call(run_meta_gen, file_path, cell_type)
    b.run(wait=False)


if __name__ == '__main__':
    main()
