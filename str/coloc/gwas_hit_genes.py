#!/usr/bin/env python3

"""
This script outputs a list of genes that have at least one SNP with pval <5e-8 in the GWAS catalog.

analysis-runner --dataset "bioheart" \
    --description "identify genomewide sig genes" \
    --access-level "test" \
    --memory='8G' \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr" \
    gwas_hit_genes.py \
    --snp-gwas-file=gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90027158.h_parsed.tsv \
    --pheno-output-name="alzheimer_GCST90027158"

"""
import click
import pandas as pd

@click.option(
    '--egenes-file',
    help='Path to the eGenes file with FINEMAP and SUSIE probabilities',
    default='gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/all_cell_types_all_genes_sig_only.tsv',
)

@click.option(
    '--snp-gwas-file',
    help='Path to the SNP GWAS file',
    default='gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv',
)

@click.option('--pheno-output-name', help='Phenotype output name', default='covid_GCST011071')
@click.command()
def main(egenes_file, snp_gwas_file, pheno_output_name):
    gwas_sig_genes = []
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
        egenes_file,
        sep='\t',
        usecols=['chr', 'pos', 'pval_meta', 'motif', 'susie_pip', 'gene', 'finemap_prob', 'celltype', 'ref_len'],
    )

    result_df_cfm = egenes
    result_df_cfm['variant_type'] = result_df_cfm['motif'].str.contains('-').map({True: 'SNV', False: 'STR'})
    result_df_cfm_str = result_df_cfm[result_df_cfm['variant_type'] == 'STR']  # filter for STRs
    result_df_cfm_str = result_df_cfm_str[result_df_cfm_str['pval_meta'] < 5e-8]  # filter for STRs with p-value < 5e-8
    result_df_cfm_str = result_df_cfm_str.drop_duplicates(
        subset=['gene'],
    )  # drop duplicates (ie pull out the distinct genes)
    result_df_cfm_str['gene'] = result_df_cfm_str['gene'].str.replace(
        '.tsv',
        '',
        regex=False,
    )  # remove .tsv from gene names (artefact of the data file)

    for gene in result_df_cfm_str['gene']:
        chrom = result_df_cfm_str[result_df_cfm_str['gene'] == gene]['chr'].iloc[0]

        # extract the coordinates for the cis-window (gene +/- 100kB)
        gene_table = var_table[var_table['gene_ids'] == gene]
        start = float(gene_table['start'].astype(float)) - 100000
        end = float(gene_table['end'].astype(float)) + 100000
        chrom = gene_table['chr'].iloc[0]
        hg38_map_chr = hg38_map[hg38_map['chromosome'] == (chrom)]
        hg38_map_chr_start = hg38_map_chr[hg38_map_chr['position'] >= start]
        hg38_map_chr_start_end = hg38_map_chr_start[hg38_map_chr_start['position'] <= end]
        if hg38_map_chr_start_end.empty:
            print('No SNP GWAS data for ' + gene + ' in the cis-window: skipping....')
            continue
        # check if the p-value column contains at least one value which is <=5e-8:
        if hg38_map_chr_start_end['p_value'].min() > 5e-8:
            print('No significant SNP GWAS data for ' + gene + ' in the cis-window: skipping....')
            continue
        gwas_sig_genes.append(gene)

    # write list to a csv file
    output_file = 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/' + pheno_output_name + '_gwas_sig_genes.csv'
    output_df = pd.DataFrame({'gene': gwas_sig_genes})
    output_df.to_csv(output_file, index=False)
    print('Output file saved:', output_file)


if __name__ == '__main__':
    main()
