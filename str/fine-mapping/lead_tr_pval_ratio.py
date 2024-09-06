#!/usr/bin/env python3

"""

This script calculates
1) ratio of p-vals between the lead TR and the lead SNV in the cis window.
2) pairwise LD (Pearson correlation) between the lead TR  and the lead SNV in the cis window.


Workflow:
2) Extract the genotypes for the lead eSTR and lead SNV in the cis window .
3) Calculate pairwise correlation of the lead eSTR and the lead SNV.
4) Save results of the top correlated variant in a TSV file.

analysis-runner --dataset "bioheart" \
    --description "Calculate LD for fine-mapped eSTRs with GWAS variants" \
    --access-level "test" \
    --memory='8G' \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr/coloc-ld/fm_strs_only/v4" \
    lead_tr_pval_ratio.py \
    --fm-csv=gs://cpg-bioheart-test/str/associatr/fine_mapping/estrs_lead_filtered.csv


"""


import click
import pandas as pd

from cpg_utils.hail_batch import get_batch
from cpg_utils import to_path


def ld_parser(
    fm_chrom,
    snp_input,
    str_input,
):
    import numpy as np
    from cyvcf2 import VCF

    max_corr_master_df = pd.DataFrame()

    for index,row in fm_chrom.iterrows():
        chrom = row['chr']
        pos = row['pos']
        end = row['end']
        motif = row['motif_x']
        cell_type = row['cell_type']
        gene = row['gene_name']
        lead_str_locus = f'{chrom}_{pos}_{motif}'
        lead_str_coord = f'{chrom}:{pos}-{end}'
        lead_str_pval = row['pval_meta']

        # get the lead SNV for that gene
        meta_results = pd.read_csv(f'gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results/{cell_type}/{chrom}/{gene}_100000bp_meta_results.tsv', sep = '\t')
        snv_meta_results = meta_results[meta_results['motif'].str.contains('-')] # filter for SNVs
        lead_snv = snv_meta_results[snv_meta_results['pval_meta'] == snv_meta_results['pval_meta'].min()] # get the lead SNV
        lead_snv_coord = chrom + ':' + str(lead_snv.iloc[0]['pos'])
        lead_snv_motif = lead_snv.iloc[0]['motif']


        pval_ratio = lead_str_pval/lead_snv.iloc[0]['pval_meta']


        # Extract the genotypes for SNVs in the cis window
        snp_df = pd.DataFrame(columns=['individual'])
        snp_vcf = VCF(snp_input['vcf'])
        snp_df['individual'] = snp_vcf.samples

        print('Starting to subset SNP VCF for window...')
        for variant in snp_vcf(lead_snv_coord):
            if str(variant.INFO.get('RU')) == lead_snv_motif:
                gt = variant.gt_types  # extracts GTs as a numpy array
                gt[gt == 3] = 2
                snp = variant.CHROM + '_' + str(variant.POS) + '_' + lead_snv_motif
                df_to_append = pd.DataFrame(gt, columns=[snp])  # creates a temp df to store the GTs for one locus
                snp_df = pd.concat([snp_df, df_to_append], axis=1)
                break
        print('Finished reading SNP VCF')

        # Extract the genotypes for STRs in the cis window
        str_df = pd.DataFrame(columns=['individual'])
        str_vcf = VCF(str_input['vcf'])
        str_df['individual'] = str_vcf.samples
        print('Starting to subset STR VCF for window...')
        for variant in str_vcf(lead_str_coord):
            if str(variant.INFO.get('RU')) == motif:

                genotypes = variant.format('REPCN')
                # Replace '.' with '-99/-99' to handle missing values
                genotypes = np.where(genotypes == '.', '-99/-99', genotypes)

                # Split each element by '/'
                split_genotypes = [genotype.split('/') for genotype in genotypes]

                # Convert split_genotypes into a numpy array for easier manipulation
                split_genotypes_array = np.array(split_genotypes)

                # Convert the strings to integers and sum them row-wise
                sums = np.sum(split_genotypes_array.astype(int), axis=1)
                # set dummy -198 value to np.nan
                sums = np.where(sums == -198, np.nan, sums)
                snp = variant.CHROM + '_' + str(variant.POS) + '_' + str(variant.INFO.get('RU'))

                str_df[snp] = sums
                break

        # merge the STR and SNP GT dfs together
        merged_df = pd.merge(str_df, snp_df, on='individual')

        # merged_df has only two columns - calculate the correlation
        correlation = merged_df.drop(columns = 'individual').corr()
        correlation.to_csv(
        f'gs://cpg-bioheart-test-analysis/str/associatr/fm_strs/pval_ratio/{cell_type}/{chrom}/{cell_type}_{chrom}_corr_table.tsv',
        sep='\t',
        index=False,
    )

        # save correlation and pval ratio to a df
        results_df = pd.DataFrame({
        'correlation': [correlation],
        'p_value': [pval_ratio],
        'lead_str_locus': [lead_str_locus],
        'lead_snv_locus': [lead_snv_coord],
        'cell_type': [cell_type],
            })


        max_corr_master_df = pd.concat([max_corr_master_df, results_df], axis=0)

    max_corr_master_df.to_csv(
        f'gs://cpg-bioheart-test-analysis/str/associatr/fm_strs/pval_ratio/{cell_type}/{chrom}/{cell_type}_{chrom}_corr.tsv',
        sep='\t',
        index=False,
    )


@click.option('--snp-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/tob_freeze_1/bgzip_tabix/v4')
@click.option('--str-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1-chr-specific')
@click.option('--fm-csv', required=True, help='Fine-mapped eSTRs CSV file path')
@click.command()
def main(fm_csv, snp_vcf_dir, str_vcf_dir):
    b = get_batch(name='Calculate LD for fine-mapped eSTRs with GWAS variants')
    fm = pd.read_csv(fm_csv)
    fm = fm.drop_duplicates(subset=['chr', 'pos', 'end', 'motif_x', 'cell_type'])
    for cell_type in fm['cell_type'].unique():
        fm_cell_type = fm[fm['cell_type'] == cell_type]
        for chrom in fm_cell_type['chr'].unique():
            if to_path(f'gs://cpg-bioheart-test-analysis/str/associatr/fm_strs/pval_ratio/{cell_type}/{chrom}/{cell_type}_{chrom}_corr.tsv').exists():
                print(f'File already exists for {cell_type} and {chrom}')
                continue

            ld_job = b.new_python_job(
                f'LD calc for {chrom}; {cell_type}',
            )
            ld_job.cpu(4)
            ld_job.storage('10G')


            fm_cell_type_chrom = fm_cell_type[fm_cell_type['chr'] == chrom]

            snp_vcf_path = f'{snp_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'
            str_vcf_path = f'{str_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'

            snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_path, 'tbi': snp_vcf_path + '.tbi'})
            str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'})

            ld_job.call(
                ld_parser,
                fm_cell_type_chrom,
                snp_input,
                str_input,

            )
            break # testing only
        break # testing only

    b.run(wait=False)


if __name__ == '__main__':
    main()