#!/usr/bin/env python3

"""
This script runs R's meta package to generate pooled effect sizes for each eQTL.
    --image australia-southeast1-docker.pkg.dev/cpg-common/images-dev/r-meta:v1 \
analysis-runner --dataset "bioheart" --description "meta results runner" --access-level "test" \
    --output-dir "str/associatr/tob_n1055_and_bioheart_n990" \
    meta_runner.py --results-dir-1=gs://cpg-bioheart-test/str/associatr/tob_n1055/results/v1 \
    --results-dir-2=gs://cpg-bioheart-test/str/associatr/bioheart_n990/results/v1 \
    --gene-list-dir-1=gs://cpg-bioheart-test/str/associatr/tob_n1055/input_files/scRNA_gene_lists/1_min_pct_cells_expressed \
    --gene-list-dir-2=gs://cpg-bioheart-test/str/associatr/bioheart_n990/input_files/scRNA_gene_lists/1_min_pct_cells_expressed \
    --cell-types=CD8_TEM --chromosomes=chr1
"""
import json

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def run_meta_gen(input_dir_1, input_dir_2, cell_type, chr, gene):
    """
    Run meta-analysis using R's meta package
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(meta)')
    ro.r('library(tidyverse)')
    d1 = pd.read_csv(f'{input_dir_1}/{cell_type}/{chr}/{gene}_100000bp.tsv', sep='\t')
    d2 = pd.read_csv(f'{input_dir_2}/{cell_type}/{chr}/{gene}_100000bp.tsv', sep='\t')

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
        f'r_df1 = r_df1 %>% select(chrom, pos, coeff_{cell_type}_{chr}_{gene},se_{cell_type}_{chr}_{gene}, motif, period, ref_len)',
    )
    ro.r(
        f'r_df2 = r_df2 %>% select(chrom, pos, coeff_{cell_type}_{chr}_{gene},se_{cell_type}_{chr}_{gene}, motif, period, ref_len)',
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
    coeff_meta = numeric(),
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
    m.gen = metagen(result_df$coeff, result_df$se, method.tau = "REML", random = TRUE,method.random.ci = "HK")


    new_entry = data.frame(
        chr = df[i, "chrom"],
        pos = df[i, "pos"],
        coeff_meta = m.gen$TE.random,
        se_meta = m.gen$seTE.random,
        pval_meta = m.gen$pval.random,
        lowerCI_meta = m.gen$lower.random,
        upperCI_meta = m.gen$upper.random
    )

    meta_df = rbind(meta_df, new_entry)
    }
    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_meta_df = ro.conversion.get_conversion().rpy2py(ro.r('meta_df'))

    pd_meta_df.to_csv(
        f'{output_path(f"meta_results/{cell_type}/{chr}/{gene}_100000bp_meta_results.tsv")}',
        sep='\t',
        index=False,
    )


@click.option('--results-dir-1', help='GCS path directory to the raw associatr results for cohort 1')
@click.option('--results-dir-2', help='GCS path directory to the raw associatr results for cohort 2')
@click.option('--gene-list-dir-1', help='GCS path directory to the gene list for cohort 1')
@click.option('--gene-list-dir-2', help='GCS path directory to the gene list for cohort 2')
@click.option('--cell-types', help='cell type')
@click.option('--chromosomes', help='chromosomes')
@click.command()
def main(results_dir_1, results_dir_2, gene_list_dir_1, gene_list_dir_2, cell_types, chromosomes):
    """
    Compute meta-analysis using R's meta package
    """

    for cell_type in cell_types.split(','):
        for chromosome in chromosomes.split(','):
            # get the intersection of genes tested in both cohorts
            gene_file_path_1 = to_path(f'{gene_list_dir_1}/{cell_type}/{chromosome}_{cell_type}_gene_list.json')
            gene_file_path_2 = to_path(f'{gene_list_dir_2}/{cell_type}/{chromosome}_{cell_type}_gene_list.json')
            with open(gene_file_path_1) as f:
                genes_1 = json.load(f)
            with open(gene_file_path_2) as g:
                genes_2 = json.load(g)
            intersected_genes =  list(set(genes_1) & set(genes_2))
            for gene in intersected_genes:
                j = get_batch(name=f'compute_meta_{cell_type}_{chromosome}_{gene}').new_python_job()
                j.cpu(1)
                j.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/r-meta:v1')
                j.call(run_meta_gen, results_dir_1, results_dir_2, cell_type, chromosome, gene)
    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
