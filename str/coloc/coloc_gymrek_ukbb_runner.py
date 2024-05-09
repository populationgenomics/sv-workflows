#!/usr/bin/env python3

"""
This script performs colocalisation analysis betweeen eGenes identified by pseudobulk STR analysis and GWAS signals from UKBB Gymrek STR catalogs.
1) Identify eGenes (FDR <5%) from the STR analysis
2) Extract the STR GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene (if STR eQTL data is available)
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by STR analysis" \
    --access-level "test" \
    --output-dir "str/associatr" \
    coloc_gymrek_ukbb_runner.py \
    --celltypes "CD4_TCM"
"""
import gzip

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def coloc_runner(gwas, eqtl_file_path, celltype, pheno):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(coloc)')
    ro.r('library(tidyverse)')

    with (ro.default_converter + pandas2ri.converter).context():
        gwas_r = ro.conversion.get_conversion().py2rpy(gwas)
    print('loaded in gwas_r')
    ro.globalenv['gwas_r'] = gwas_r
    ro.r(
        '''
    print(names(gwas_r))
    gwas_r$pvalues = gwas_r$p_value
    gwas_r$varbeta = (gwas_r$standard_error)**2
    gwas_r$position = gwas_r$start_pos..hg38. -1
    gwas_r$snp = paste('s', gwas_r$position, sep = '')

    gwas_r = gwas_r %>% select(beta, varbeta, position,snp)
    gwas_r = gwas_r%>% as.list()
    gwas_r$type = 'quant'
    gwas_r$sdY = 1

     ''',
    )
    eqtl = pd.read_csv(eqtl_file_path, sep='\t')
    gene = eqtl_file_path.split('/')[-1].split('_')[0]
    with (ro.default_converter + pandas2ri.converter).context():
        eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
    print('loaded in eqtl_r')
    ro.globalenv['eqtl_r'] = eqtl_r
    ro.globalenv['gene'] = gene
    ro.r(
        '''
    eqtl_r$beta = eqtl_r$coeff_meta
    eqtl_r$varbeta = eqtl_r$se_meta**2
    eqtl_r$position = eqtl_r$pos
    eqtl_r$snp = paste('s', eqtl_r$position, sep = '')
    eqtl_r = eqtl_r%>% select(beta, varbeta, position, snp)
    eqtl_r = eqtl_r %>% distinct(snp, .keep_all = TRUE)


    eqtl_r = eqtl_r %>% as.list()
    eqtl_r$type = 'quant'
    eqtl_r$sdY = 1
    my.res <- coloc.abf(dataset1=gwas_r,
                    dataset2=eqtl_r)
    p_df <- data.frame(gene, my.res$summary[6])
    names(p_df) <- c('gene', 'PP.H4.abf')
    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_p4_df = ro.conversion.get_conversion().rpy2py(ro.r('p_df'))
    print('converted back to pandas df')

    # add cell type annotation to df
    pd_p4_df['celltype'] = celltype

    # write to GCS
    pd_p4_df.to_csv(
        f'{output_path(f"coloc/gymrek-ukbb-{pheno}/{celltype}/{gene}_100kb.tsv")}',
        sep='\t',
        index=False,
    )


@click.option(
    '--egenes-dir',
    help='Path to the eGenes dir',
    default='gs://cpg-bioheart-test-analysis/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat',
)
@click.option(
    '--str-cis-dir',
    help='Path to the directory containing the STR eQTL cis results',
    default='gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results',
)
@click.option('--celltypes', help='Cell type for which the eGenes were identified', default='CD4_TCM')
@click.option(
    '--var-annotation-file',
    help='Gene annotation file path',
    default='gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv',
)
@click.option('--pheno', help='Phenotype to use for coloc', default='alanine_aminotransferase')
@click.command()
def main(str_cis_dir, egenes_dir, celltypes, var_annotation_file, pheno):
    # read in gene annotation file
    var_table = pd.read_csv(var_annotation_file)

    for phenotype in pheno.split(','):
        gwas_file = to_path(f'gs://cpg-bioheart-test/str/gymrek-ukbb-str-gwas-catalogs/gymrek-ukbb-str-gwas-catalogs/white_british_{phenotype}_str_gwas_results.tab.gz')
        with gzip.open(gwas_file, 'rb') as f:
            gwas = pd.read_csv(
                f,
                sep='\t',
                usecols=['chromosome', 'beta', 'standard_error', 'p_value', 'repeat_unit', 'start_pos (hg38)'],
            )

        for celltype in celltypes.split(','):
            egenes_file_path = f'{egenes_dir}/{celltype}_qval.tsv'
            # read in eGenes file
            egenes = pd.read_csv(egenes_file_path, sep='\t')
            egenes = egenes[egenes['qval'] < 0.05]  # filter for eGenes with FDR<5%

            #for gene in egenes['gene_name']:
            for gene in ['ENSG00000107771']:
                # extract the coordinates for the cis-window (gene +/- 100kB)
                gene_table = var_table[var_table['gene_ids'] == gene]
                start = float(gene_table['start'].astype(float)) - 100000
                end = float(gene_table['end'].astype(float)) + 100000
                chrom = gene_table['chr'].iloc[0][3:]
                hg38_map_chr = gwas[gwas['chromosome'] == int(chrom)]
                hg38_map_chr_start = hg38_map_chr[hg38_map_chr['start_pos (hg38)'] >= start]
                hg38_map_chr_start_end = hg38_map_chr_start[hg38_map_chr_start['start_pos (hg38)'] <= end]

                if hg38_map_chr_start_end.empty:
                    print('No STR GWAS data for ' + gene + ' in the cis-window: skipping....')
                    continue
                print('Extracted SNP GWAS data for ' + gene)

                # run coloc
                b = get_batch()
                coloc_job = b.new_python_job(
                    f'Coloc for {gene}: {celltype}',
                )
                coloc_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/r-meta:2.0')
                coloc_job.call(
                    coloc_runner,
                    hg38_map_chr_start_end,
                    f'{str_cis_dir}/{celltype}/chr{chrom}/{gene}_100000bp_meta_results.tsv',
                    celltype,
                    phenotype,
                )

    b.run(wait=False)


if __name__ == '__main__':
    main()
