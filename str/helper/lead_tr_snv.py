#!/usr/bin/env python3

"""
Calculates LD between lead variant (TR/SNV) and next significant SNV (proxy SNV) for each gene x cell type combination.

 analysis-runner  --dataset "tenk10k" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr/final_freeze/meta_fixed" \
python3 lead_tr_snv.py

"""
import json

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch
from cpg_utils.hail_batch import output_path


def genes_parser(
    chromosome,
    cell_type,
    snp_input,
    str_input,
):
    """ """
    import numpy as np
    from cyvcf2 import VCF

    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    max_corr_master_df = pd.DataFrame()

    chromosome = f'chr{chromosome}'
    genes = list(to_path(f'gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/trs_snps/rm_str_indels_dup_strs_v3/{cell_type}/{chromosome}').rglob('*.tsv'))
    for gene_file in genes:
        eqtl_results = pd.read_csv(
            gene_file,
            sep='\t',)
        gene_name = str(gene_file).split('/')[-1].split('_')[0]
        # get row(s) with minimum p-value
        min_pval = eqtl_results['pval_meta_fixed'].min()
        smallest_pval_rows = eqtl_results[eqtl_results['pval_meta_fixed'] == min_pval]
        # check if all rows are SNPs:
        at_least_one_lead_SNV = smallest_pval_rows['motif'].str.contains('-').any()
        if at_least_one_lead_SNV:
            snv_meta_results = eqtl_results[eqtl_results['motif'].str.contains('-')]  # filter for SNVs
            lead_snv = snv_meta_results[
                snv_meta_results['pval_meta_fixed'] == snv_meta_results['pval_meta_fixed'].min()
            ]  # get the lead SNV
            lead_snv_coord = chromosome + ':' + str(lead_snv.iloc[0]['pos'])
            lead_snv_motif = lead_snv.iloc[0]['motif']

            proxy_snv_results = snv_meta_results[snv_meta_results.index != lead_snv.index[0]]
            proxy_snv = proxy_snv_results[proxy_snv_results['pval_meta_fixed'] == proxy_snv_results['pval_meta_fixed'].min()]
            proxy_snv_coord = f"{chromosome}:{proxy_snv.iloc[0]['pos']}"
            proxy_snv_motif = proxy_snv.iloc[0]['motif']
            distance = abs(int(lead_snv.iloc[0]['pos']) - int(proxy_snv.iloc[0]['pos']))

            # Extract the genotypes for lead SNV in the cis window
            lead_df = pd.DataFrame(columns=['individual'])
            snp_vcf_1 = VCF(snp_input['vcf'])
            lead_df['individual'] = snp_vcf_1.samples

            for variant in snp_vcf_1(lead_snv_coord):
                if str(variant.INFO.get('RU')) == lead_snv_motif:
                    gt = variant.gt_types  # extracts GTs as a numpy array
                    gt[gt == 3] = 2
                    snp = variant.CHROM + '_' + str(variant.POS) + '_' + lead_snv_motif
                    df_to_append = pd.DataFrame(gt, columns=[snp])  # creates a temp df to store the GTs for one locus
                    lead_df = pd.concat([lead_df, df_to_append], axis=1)
                    break
                print('Finished reading SNP VCF')

        else:
            lead_tr = eqtl_results[eqtl_results['pval_meta_fixed'] == min_pval]
            lead_tr_coord = chromosome + ':' + str(lead_tr.iloc[0]['pos'])
            lead_tr_motif = lead_tr.iloc[0]['motif']

            proxy_snv_results = eqtl_results[eqtl_results.index != lead_tr.index[0]]  # filter out lead tr
            proxy_snv_results = proxy_snv_results[proxy_snv_results['motif'].str.contains('-')]  # filter for SNVs
            proxy_snv = proxy_snv_results[proxy_snv_results['pval_meta_fixed'] == proxy_snv_results['pval_meta_fixed'].min()]
            proxy_snv_coord = f"{chromosome}:{proxy_snv.iloc[0]['pos']}"
            proxy_snv_motif = proxy_snv.iloc[0]['motif']
            distance = abs(int(lead_tr.iloc[0]['pos']) - int(proxy_snv.iloc[0]['pos']))

            # Extract the genotypes for STRs in the cis window
            lead_df = pd.DataFrame(columns=['individual'])
            str_vcf = VCF(str_input['vcf'])
            lead_df['individual'] = str_vcf.samples
            print('Starting to subset STR VCF for window...')
            for variant in str_vcf(lead_tr_coord):
                if str(variant.INFO.get('RU')) == lead_tr_motif:
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

                    lead_df[snp] = sums
                    break

        # Extract the genotypes for proxy SNV in the cis window
        proxy_df = pd.DataFrame(columns=['individual'])
        snp_vcf_2 = VCF(snp_input['vcf'])
        proxy_df['individual'] = snp_vcf_2.samples

        for variant in snp_vcf_2(proxy_snv_coord):
            if str(variant.INFO.get('RU')) == proxy_snv_motif:
                gt = variant.gt_types
                gt[gt == 3] = 2
                snp = variant.CHROM + '_' + str(variant.POS) + '_' + proxy_snv_motif
                df_to_append = pd.DataFrame(gt, columns=[snp])
                proxy_df = pd.concat([proxy_df, df_to_append], axis=1)
                break
            print('Finished reading SNP VCF')

        # calculate LD
        # merge the STR and SNP GT dfs together
        merged_df = lead_df.merge(proxy_df, on='individual')

        # merged_df has only two columns - calculate the correlation
        try:
            correlation = merged_df.drop(columns='individual').corr().iloc[0, 1]
        except:
            print(correlation)
            print(gene_name)
            continue


        # save correlation and pval ratio to a df
        results_df = pd.DataFrame(
            {
                'correlation': [correlation],
                'proxy_pval': [proxy_snv.iloc[0]['pval_meta_fixed']],
                'lead_pval': [min_pval],
                'lead_snv_boolean': [at_least_one_lead_SNV],
                'cell_type': [cell_type],
                'gene': [gene_name],
                'distance': [distance],
            },
        )

        max_corr_master_df = pd.concat([max_corr_master_df, results_df], axis=0)
        max_corr_master_df.to_csv(
            output_path(
                f'lead_tr_snv_proxy/{cell_type}/{chromosome}/{cell_type}_{chromosome}_summ_stats.tsv',
                'analysis',
            ),
            sep='\t',
            index=False,
        )


@click.option('--snp-vcf-dir', default='gs://cpg-bioheart-test/tenk10k/str/associatr/common_variant_snps')
@click.option('--str-vcf-dir', default='gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific')
@click.command()
def main(snp_vcf_dir, str_vcf_dir):
    """
    Calculates LD between lead variant (TR/SNV) and next significant SNV (proxy SNV) for each gene x cell type combination.
    """
    b = get_batch(name='Lead TR, SNV, and proxy SNV data extraction')

    celltypes = 'gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive'
    celltypes = celltypes.split(',')
    for cell_type in celltypes:
        for chrom in range(1, 23):
            if to_path(
                output_path(f'lead_tr_snv_proxy/{cell_type}/{chrom}/{cell_type}_{chrom}_summ_stats.tsv',
            )).exists():
                print(f'File already exists for {cell_type} and {chrom}')
                continue
            j = b.new_python_job(
                name=f'Get pvals/LD of lead variant and closest SNP proxy {cell_type}: {chrom}',
            )
            snp_vcf_path = f'{snp_vcf_dir}/hail_filtered_chr{chrom}.vcf.bgz'
            str_vcf_path = f'{str_vcf_dir}/hail_filtered_chr{chrom}.vcf.bgz'

            snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_path, 'tbi': snp_vcf_path + '.tbi'})
            str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'})
            j.call(
                genes_parser,
                chrom,
                cell_type,
                snp_input,
                str_input,
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments
