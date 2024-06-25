#!/usr/bin/env python3

"""
This script performs SNP-only colocalisation analysis betweeen eGenes identified by pseudobulk STR analysis and GWAS signals.
Assumes that the SNP GWAS data has been pre-processed with the following columns: 'chromosome', 'position' (hg38 bp), 'snp'(chromosome_position_refallele_effectallele), 'beta', 'varbeta'

1) Identify eGenes where at least one STR has a fine-map causal probability of at least 80% by both SUSIE and FINEMAP.
2) Extract the SNP GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by STR analysis" \
    --access-level "test" \
    --output-dir "str/associatr" \
    coloc_runner.py
    --celltypes "CD4_TCM"

"""

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, image_path


def coloc_runner(gwas, eqtl_file_path, celltype, pheno_output_name):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(coloc)')
    ro.r('library(tidyverse)')

    with (ro.default_converter + pandas2ri.converter).context():
        gwas_r = ro.conversion.get_conversion().py2rpy(gwas)
    ro.globalenv['gwas_r'] = gwas_r
    ro.r(
        '''
    gwas_r = gwas_r %>% select(beta, varbeta, position,snp)
    gwas_r = gwas_r %>% distinct(snp, .keep_all = TRUE)
    gwas_r$type = 'cc'
    gwas_r = gwas_r%>% as.list()
    ''',
    )
    eqtl = pd.read_csv(
        eqtl_file_path,
        sep='\t',
        names=[
            'chrom',
            'pos',
            'alleles',
            'n_samples_tested',
            'locus_filtered',
            'p-val',
            'beta',
            'se',
            'r2',
            'motif',
            'period',
            'ref-len',
            'allele-freq',
        ],
        header=0,
    )
    eqtl['snp'] = eqtl['chrom'] + '_' + eqtl['pos'].astype(str) + '_' + eqtl['motif']
    eqtl['snp'] = eqtl['snp'].str.replace('-', '_', regex=False)
    gene = eqtl_file_path.split('/')[-1].split('_')[2]
    with (ro.default_converter + pandas2ri.converter).context():
        eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
    ro.globalenv['eqtl_r'] = eqtl_r
    ro.globalenv['gene'] = gene
    ro.r(
        '''
    eqtl_r = eqtl_r %>% filter(!is.na(beta))
    eqtl_r = eqtl_r %>% distinct(snp, .keep_all = TRUE)
    eqtl_r$varbeta = eqtl_r$se**2
    eqtl_r$position = eqtl_r$pos
    eqtl_r = eqtl_r %>% select(beta, varbeta, position, snp)

    eqtl_r = eqtl_r %>% as.list()
    eqtl_r$type = 'quant'
    eqtl_r$sdY = 1


    my.res <- coloc.abf(dataset1=gwas_r,
                    dataset2=eqtl_r)

    p_df <- data.frame(gene, my.res$summary[1:6])
    names(p_df) <- c('gene', 'nsnps_coloc_tested','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
    ''',
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_p4_df = ro.conversion.get_conversion().rpy2py(ro.r('p_df'))

    # write to GCS
    pd_p4_df.to_csv(
        output_path(f"coloc-snp-only/{pheno_output_name}/{celltype}/{gene}_100kb.tsv", 'analysis'),
        sep='\t',
        index=False,
    )


@click.option(
    '--egenes',
    help='Path to the eGenes file with FINEMAP and SUSIE probabilities',
    default='gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/all_cell_types_all_genes_sig_only.tsv',
)
@click.option(
    '--snp-cis-dir',
    help='Path to the directory containing the SNP cis results',
    default='gs://cpg-bioheart-test/str/associatr/common_variants_snps/tob_n1055_and_bioheart_n990/meta_results/meta_results',
)
@click.option(
    '--snp-gwas-file',
    help='Path to the SNP GWAS file',
    default='gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv',
)
@click.option('--celltypes', help='Cell types to run', default='ASDC')
@click.option('--pheno-output-name', help='Phenotype output name', default='covid_GCST011071')
@click.command()
def main(snp_cis_dir, egenes, celltypes, snp_gwas_file, pheno_output_name):
    # read in gene annotation file
    var_table = pd.read_csv(
        'gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv',
    )
    hg38_map = pd.read_csv(
        snp_gwas_file,
        sep='\t',
    )

    # read in eGenes file
    egenes = pd.read_csv(
        egenes,
        sep='\t',
        usecols=['chr', 'pos', 'pval_meta', 'motif', 'susie_pip', 'gene', 'finemap_prob', 'celltype', 'ref_len'],
    )
    result_df_80 = egenes[
        (egenes['susie_pip'] >= 0.8) & (egenes['finemap_prob'] >= 0.8)
    ]  # filter for variants where the causal probability is at least 80% by both SUSIE and FINEMAP
    result_df_cfm = result_df_80[result_df_80['pval_meta'] < 1e-10]  # filter for variants with p-value < 1e-10
    result_df_cfm['variant_type'] = result_df_cfm['motif'].str.contains('-').map({True: 'SNV', False: 'STR'})
    result_df_cfm_str = result_df_cfm[result_df_cfm['variant_type'] == 'STR']  # filter for STRs
    result_df_cfm_str = result_df_cfm_str.drop_duplicates(
        subset=['gene', 'celltype'],
    )  # drop duplicates (ie pull out the distinct genes in each celltype)
    result_df_cfm_str['gene'] = result_df_cfm_str['gene'].str.replace(
        '.tsv', '', regex=False,
    )  # remove .tsv from gene names (artefact of the data file)
    b = get_batch(name='Run coloc')

    for celltype in celltypes.split(','):
        result_df_cfm_str_celltype = result_df_cfm_str[
            result_df_cfm_str['celltype'] == celltype
        ]  # filter for the celltype of interest
        for gene in result_df_cfm_str_celltype['gene']:
            chr = result_df_cfm_str_celltype[result_df_cfm_str_celltype['gene'] == gene]['chr'].iloc[0]
            if to_path(f'{snp_cis_dir}/{celltype}/{chr}/{gene}_100000bp_meta_results.tsv').exists():
                print('Cis results for ' + gene + ' exist: proceed with coloc')

                # extract the coordinates for the cis-window (gene +/- 100kB)
                gene_table = var_table[var_table['gene_ids'] == gene]
                start = float(gene_table['start'].astype(float)) - 100000
                end = float(gene_table['end'].astype(float)) + 100000
                chr = gene_table['chr'].iloc[0]
                hg38_map_chr = hg38_map[hg38_map['chromosome'] == (chr)]
                hg38_map_chr_start = hg38_map_chr[hg38_map_chr['position'] >= start]
                hg38_map_chr_start_end = hg38_map_chr_start[hg38_map_chr_start['position'] <= end]
                if hg38_map_chr_start_end.empty:
                    print('No SNP GWAS data for ' + gene + ' in the cis-window: skipping....')
                    continue
                print('Extracted SNP GWAS data for ' + gene)

                # run coloc
                coloc_job = b.new_python_job(
                    f'Coloc for {gene}: {celltype}',
                )
                f'{snp_cis_dir}/{celltype}/{chr}/{gene}_100000bp_meta_results.tsv'
                coloc_job.image(image_path('r-meta'))
                coloc_job.call(
                    coloc_runner,
                    hg38_map_chr_start_end,
                    f'{snp_cis_dir}/{celltype}/{chr}/{gene}_100000bp_meta_results.tsv',
                    celltype,
                )

            else:
                print('No cis results for ' + gene + ' exist: skipping....')

    b.run(wait=False)


if __name__ == '__main__':
    main()
