#!/usr/bin/env python3

"""
This script performs colocalisation analysis betweeen eGenes identified by pseudobulk STR analysis and GWAS signals.

1) Identify eGenes (FDR <5%) from the STR analysis
2) Extract the SNP GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene (if SNP eQTL data is available
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by STR analysis" \
    --access-level "test" \
    --output-dir "str/associatr" \
    coloc_runner.py --egenes "gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat/CD4_TCM_qval.tsv" \
    --snp-cis-dir "gs://cpg-bioheart-test/saige-qtl/bioheart_n990/v1" \
    --celltype "CD4_TCM"


"""

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def coloc_runner(gwas, eqtl_file_path, celltype):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(coloc)')
    ro.r('library(todyverse)')

    with (ro.default_converter + pandas2ri.converter).context():
        gwas_r = ro.conversion.get_conversion().py2rpy(gwas)
    ro.globalenv['gwas_r'] = gwas_r
    ro.r(
        '''
    gwas_r$MAF = (gwas_r$FRQ_A_12882*12882 + gwas_r$FRQ_U_21770*21770)/(12882+21770)
    gwas_r$pvalues = gwas_r$P
    gwas_r$or = ggwas_rwas$OR
    gwas_r$beta = log(gwas_r$or)
    gwas_r$varbeta = (gwas_r$SE/gwas$or)**2
    gwas_r$position = gwas_r$BP
    gwas_r$snp = paste('s', gwas_r$BP, sep = '')

    gwas_r = gwas_r %>% select(beta, varbeta, position,snp,MAF)
    gwas_r = gwas_r%>% as.list()
    gwas_r$type = 'cc'
    gwas_r$N = 34652
    gwas_r$s = 0.33

     '''
    )
    eqtl = pd.read_csv(eqtl_file_path, sep='\t')
    gene = eqtl_file_path.split('/')[-1].split('_')[1]
    with (ro.default_converter + pandas2ri.converter).context():
        eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
    ro.globalenv['eqtl_r'] = eqtl_r
    ro.globalenv['gene'] = gene
    ro.r(
        '''
    eqtl_r = eqtl_r %>% filter(!is.na(BETA))
    eqtl_r = eqtl_r %>% distinct(POS, .keep_all = TRUE)
    eqtl_r$beta = eqtl_r$BETA
    eqtl_r$varbeta = eqtl_r$SE**2
    eqtl_r$position = eqtl_r$POS
    eqtl_r$snp = paste('s', eqtl_r$position, sep = '')
    eqtl_r = eqtl_r %>% select(beta, varbeta, position, snp)

    eqtl_r = eqtl_r %>% as.list()
    eqtl_r$type = 'quant'
    eqtl_r$sdY = 1


    my.res <- coloc.abf(dataset1=gwas_r,
                    dataset2=eqtl_r)

    p_df <- data.frame(gene, my.res$summary[6])
    rownames(p_df) <- c('gene', 'PP.H4.abf')
    '''
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_p4_df = ro.conversion.get_conversion().rpy2py(ro.r('p_df'))

    # write to GCS
    pd_p4_df.to_csv(
        f'{output_path(f"coloc/{celltype}/{gene}_100kb.tsv")}',
        sep='\t',
        index=False,
    )


@click.option(
    '--egenes',
    help='Path to the eGenes file',
    default='gs://cpg-bioheart-test-analysis/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat/CD4_TCM_qval.tsv',
)
@click.option(
    '--snp-cis-dir',
    help='Path to the directory containing the SNP cis results',
    default='gs://cpg-bioheart-test-analysis/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat/snp_cis',
)
@click.option('--celltype', help='Cell type for which the eGenes were identified', default='CD4_TCM')
@click.command()
def main(snp_cis_dir, egenes, celltype):
    # read in gene annotation file
    var_table = pd.read_csv(
        'gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv'
    )
    hg38_map = pd.read_csv(
        'hg38.EUR.IBD.gwas_info03_filtered.assoc',
        sep='\t',
        header=None,
        names=['CHR', 'BP', 'END', 'FRQ_A_12882', 'FRQ_U_21770', 'OR', 'SE', 'P'],
    )

    # read in eGenes file
    egenes = pd.read_csv(egenes, sep='\t')
    egenes = egenes[egenes['qval'] < 0.05]  # filter for eGenes with FDR<5%

    for gene in egenes['gene_name']:
        if to_path(snp_cis_dir + '/' + celltype + '_' + gene + '_cis').exists():
            print('Cis results for ' + gene + ' exist: proceed with coloc')

            # extract the coordinates for the cis-window (gene +/- 100kB)
            gene_table = var_table[var_table['gene_ids'] == gene]
            start = float(gene_table['start'].astype(float)) - 100000
            end = float(gene_table['end'].astype(float)) + 100000
            chr = gene_table['chr'].iloc[0][3:]
            hg38_map_chr = hg38_map[hg38_map['CHR'] == int(chr)]
            hg38_map_chr_start = hg38_map_chr[hg38_map_chr['BP'] >= start]
            hg38_map_chr_start_end = hg38_map_chr_start[hg38_map_chr_start['BP'] <= end]
            print('Extracted SNP GWAS data for ' + gene)

            # run coloc
            b = get_batch()
            coloc_job = b.new_python_job(
                f'Coloc for {gene}: {celltype}',
            )
            coloc_job.call(
                coloc_runner, hg38_map_chr_start_end, snp_cis_dir + '/' + celltype + '_' + gene + '_cis', celltype
            )

        else:
            print('No cis results for ' + gene + ' exist: skipping....')

        b.run(wait=False)


if __name__ == '__main__':
    main()
