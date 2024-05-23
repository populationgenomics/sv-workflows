#!/usr/bin/env python3

"""
This script performs colocalisation analysis betweeen eGenes identified by pseudobulk STR analysis and GWAS signals from IBD meta-analysis.

1) Identify eGenes (FDR <5%) from the STR analysis
2) Extract the SNP GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene - only do this if STR is still the lead signal in the joint SNP+STR eQTL analysis.
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by STR analysis" \
    --access-level "test" \
    --output-dir "str/associatr" \
    coloc_ibd_runner.py --egenes-dir "gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat" \
    --snp-cis-dir "gs://cpg-bioheart-test/str/associatr/common_variants_snps/tob_n1055_and_bioheart_n990/meta_results" \
    --celltypes "CD4_TCM"


"""

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

def coloc_runner(gene, snp_cis_dir, var_annotation_file, gwas_file, celltype):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    # read in gene annotation file
    var_table = pd.read_csv(var_annotation_file)

    hg38_map = pd.read_csv(
        gwas_file,
        sep='\t',
        header=None,
        names=['CHR', 'BP', 'END', 'FRQ_A_12882', 'FRQ_U_21770', 'OR', 'SE', 'P'],
    )

    gene_table = var_table[var_table['gene_ids'] == gene]
    chr = gene_table['chr'].iloc[0][3:]
    if to_path(f'{snp_cis_dir}/{celltype}/chr{chr}/{gene}_100000bp_meta_results.tsv').exists():
        print('Cis results for ' + gene + ' exist: proceed with coloc')

        # extract the coordinates for the cis-window (gene +/- 100kB)
        start = float(gene_table['start'].astype(float)) - 100000
        end = float(gene_table['end'].astype(float)) + 100000
        hg38_map_chr = hg38_map[hg38_map['CHR'] == int(chr)]
        hg38_map_chr_start = hg38_map_chr[hg38_map_chr['BP'] >= start]
        hg38_map_chr_start_end = hg38_map_chr_start[hg38_map_chr_start['BP'] <= end]
        hg38_map_chr_start_end['locus'] = (
            hg38_map_chr_start_end['CHR'].astype(str) + ':' + hg38_map_chr_start_end['BP'].astype(str)
        )
        hg38_map_chr_start_end[['locus']].to_csv(
            output_path(f'coloc-snp-only/lead-str-signal/ibd/{celltype}/{gene}_snp_gwas_list.csv'),
            index=False,
        )
        if hg38_map_chr_start_end.empty:
            print('No SNP GWAS data for ' + gene + ' in the cis-window: skipping....')
            return
        print('Extracted SNP GWAS data for ' + gene)
    else:
        print('No cis results for ' + gene + ' exist: skipping....')
        return

    ro.r('library(coloc)')
    ro.r('library(tidyverse)')

    with (ro.default_converter + pandas2ri.converter).context():
        gwas_r = ro.conversion.get_conversion().py2rpy(hg38_map_chr_start_end)
    ro.globalenv['gwas_r'] = gwas_r
    ro.r(
        '''
    gwas_r$MAF = (gwas_r$FRQ_A_12882*12882 + gwas_r$FRQ_U_21770*21770)/(12882+21770)
    gwas_r$pvalues = gwas_r$P
    gwas_r$or = gwas_r$OR
    gwas_r$beta = log(gwas_r$or)
    gwas_r$varbeta = (gwas_r$SE/gwas_r$or)**2
    gwas_r$position = gwas_r$BP
    gwas_r$snp = paste('s', gwas_r$BP, sep = '')

    gwas_r = gwas_r %>% select(beta, varbeta, position,snp,MAF)
    gwas_r = gwas_r %>% distinct(position, .keep_all = TRUE)
    gwas_r = gwas_r%>% as.list()
    gwas_r$type = 'cc'
    gwas_r$N = 34652
    gwas_r$s = 0.33

     ''',
    )
    eqtl = pd.read_csv(f'{snp_cis_dir}/{celltype}/chr{chr}/{gene}_100000bp_meta_results.tsv', sep='\t')
    with (ro.default_converter + pandas2ri.converter).context():
        eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
    ro.globalenv['eqtl_r'] = eqtl_r
    ro.globalenv['gene'] = gene
    ro.r(
        '''
    eqtl_r = eqtl_r %>% filter(!is.na(coeff_meta))
    eqtl_r = eqtl_r %>% distinct(pos, .keep_all = TRUE)
    eqtl_r$beta = eqtl_r$coeff_meta
    eqtl_r$varbeta = eqtl_r$se_meta**2
    eqtl_r$position = eqtl_r$pos
    eqtl_r$snp = paste('s', eqtl_r$position, sep = '')
    eqtl_r = eqtl_r %>% select(beta, varbeta, position, snp)

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

    # add cell type annotation to df
    pd_p4_df['celltype'] = celltype

    # write to GCS
    pd_p4_df.to_csv(
        f'{output_path(f"coloc-snp-only/lead-str-signal/ibd/{celltype}/{gene}_100kb.tsv")}',
        sep='\t',
        index=False,
    )


@click.option(
    '--egenes-dir',
    help='Path to the eGenes dir',
    default='gs://cpg-bioheart-test-analysis/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat',
)
@click.option(
    '--snp-cis-dir',
    help='Path to the directory containing the SNP cis results',
    default='gs://cpg-bioheart-test/saige-qtl/bioheart_n990/v2/output_files/output_files',
)
@click.option('--celltypes', help='Cell type for which the eGenes were identified', default='CD4_TCM')
@click.option(
    '--var-annotation-file',
    help='Gene annotation file path',
    default='gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv',
)
@click.option(
    '--gwas-file',
    help='Path to the GWAS file',
    default='gs://cpg-bioheart-test/str/gwas_catalog/g38.EUR.IBD.gwas_info03_filtered.assoc',
)
@click.command()
def main(snp_cis_dir, egenes_dir, celltypes, var_annotation_file, gwas_file):
    b = get_batch()

    for celltype in celltypes.split(','):
        egenes_file_path = f'{egenes_dir}/{celltype}_qval.tsv'
        # read in eGenes file
        egenes = pd.read_csv(egenes_file_path, sep='\t')
        egenes = egenes[egenes['qval'] < 0.05]  # filter for eGenes with FDR<5%

        # read in significant eSTRs file (joint analysis
        estrs = pd.read_csv('gs://cpg-bioheart-test/str/associatr/eSTRs_from_joint_analysis.csv')
        #subset to cell type
        estrs = estrs[estrs['cell_type']==celltype]
        for gene in egenes['gene_name']:
                if gene in set(estrs['gene_name']):
                    # run coloc
                    coloc_job = b.new_python_job(
                        f'Coloc for {gene}: {celltype}',
                    )
                    coloc_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images/r-meta:7.0.0')
                    coloc_job.call(
                        coloc_runner,
                        gene,
                        snp_cis_dir,
                        var_annotation_file,
                        gwas_file,
                        celltype,
                    )

    b.run(wait=False)


if __name__ == '__main__':
    main()
