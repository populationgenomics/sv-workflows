#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script writes out genes for every cell type if the gene has a lead signal that is not a SNP.
We test all genes tested, not just eGenes that pass an FDR.
Used to generate one stat in the paper.

 analysis-runner  --dataset "bioheart" --access-level "test" \
 --storage="20G" --memory='8G' \
--description "get cis and numpy" --output-dir "str/associatr/estrs" \
python3 lead_tr_snv_ld_decay_metrics.py

"""
import json

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def genes_parser(
    chromosome,
    cell_type,
    snp_input,
    str_input,
    gene
):
    """ """
    import numpy as np
    from cyvcf2 import VCF

    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    max_corr_master_df = pd.DataFrame()

    chromosome = f'chr{chromosome}'


    # for gene in genes:
    for gene in [gene]:
        try:
            eqtl_results = pd.read_csv(
                f'gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results/{cell_type}/{chromosome}/{gene}_100000bp_meta_results.tsv',
                sep='\t',
            )
        except FileNotFoundError:
            print(f'No eQTL results found for {gene}... skipping')
            continue
        # get row(s) with minimum p-value
        min_pval = eqtl_results['pval_meta'].min()
        smallest_pval_rows = eqtl_results[eqtl_results['pval_meta'] == min_pval]
        # check if all rows are SNPs:
        at_least_one_lead_SNV = smallest_pval_rows['motif'].str.contains('-').all()
        if at_least_one_lead_SNV:
            snv_meta_results = eqtl_results[eqtl_results['motif'].str.contains('-')]  # filter for SNVs
            lead_snv = snv_meta_results[
                snv_meta_results['pval_meta'] == snv_meta_results['pval_meta'].min()
            ]  # get the lead SNV
            lead_snv_coord = chromosome + ':' + str(lead_snv.iloc[0]['pos'])
            lead_snv_motif = lead_snv.iloc[0]['motif']

            proxy_snv_results = snv_meta_results[snv_meta_results.index != lead_snv.index[0]]
            proxy_snv = proxy_snv_results[proxy_snv_results['pval_meta'] == proxy_snv_results['pval_meta'].min()]

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
            lead_variant_coord = lead_snv_coord

        elif not smallest_pval_rows['motif'].str.contains('-').any():
            lead_tr = eqtl_results[eqtl_results['pval_meta'] == min_pval]
            lead_tr_coord = chromosome + ':' + str(lead_tr.iloc[0]['pos'])
            lead_tr_motif = lead_tr.iloc[0]['motif']

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
            lead_variant_coord = lead_tr_coord
        else:
            print('No lead variant found')
            continue
        # Extract max R2 for lead TR and SNVs in 10kb bins +/- 500kb from lead variant
        ## Define the bins
        lead_variant_chrom = lead_variant_coord.split(':')[0]
        lead_variant_pos = int(lead_variant_coord.split(':')[1])
        print(f'Lead variant position: {lead_variant_pos}')
        one_mb_window_start = max(lead_variant_pos - 500000, 0)
        one_mb_window_end = lead_variant_pos + 500000
        bins = [
            f'{lead_variant_chrom}:{start+1}-{start+10000}'
            for start in range(one_mb_window_start, one_mb_window_end, 10000)
        ]
        ## Iterate through the bins
        for i, bin in enumerate(bins):
            ## Extract the genotypes for SNVs in the bin
            snp_df = pd.DataFrame(columns=['individual'])
            snp_vcf = VCF(snp_input['vcf'])
            snp_df['individual'] = snp_vcf.samples
            print(bin)
            for variant in snp_vcf(bin):
                gt = variant.gt_types  # extracts GTs as a numpy array
                gt[gt == 3] = 2
                snp = variant.CHROM + '_' + str(variant.POS) + '_' + str(variant.INFO.get('RU'))
                df_to_append = pd.DataFrame(gt, columns=[snp])  # creates a temp df to store the GTs for one locus
                snp_df = pd.concat([snp_df, df_to_append], axis=1)

            # Calculate the max R2 for the lead TR and SNVs in the 10kb bins
            merged_df = lead_df.merge(snp_df, on='individual')
            # Get correlation matrix
            corr_matrix = merged_df.drop(columns='individual').corr()

            # Get maximum absolute correlation off diagonal
            try:
                # get the max absolute correlation in the first column (corresponds to lead variant) but it cant be on the diagnonal
                first_column = corr_matrix.iloc[:, 0] if hasattr(corr_matrix, 'iloc') else corr_matrix[:, 0]  # Extract the first column
                max_abs_corr = np.max(np.abs(first_column[1:]))  # Exclude the diagonal by slicing (start from index 1)
            except ValueError: # if there are no off-diagonal elements
                print(corr_matrix)
                print(f'No off-diagonal elements for {bin}')
                max_abs_corr = np.nan


            # save results for input into results_df
            # Store max correlation for this bin
            bin_midpoint = (int(bin.split(':')[1].split('-')[1]) + int(bin.split(':')[1].split('-')[0]))/2
            distance = bin_midpoint - lead_variant_pos

            # Create results DataFrame with common attributes and bin correlations
            results_df = pd.DataFrame(
                {
                    'lead_pval': [min_pval],
                    'lead_snv_boolean': [at_least_one_lead_SNV],
                    'cell_type': [cell_type],
                    'gene': [gene],
                    'distance': [distance],
                    'max_abs_corr': [max_abs_corr],
                },
            )

            # Append results to master DataFrame
            max_corr_master_df = pd.concat([max_corr_master_df, results_df], axis=0)
    max_corr_master_df.to_csv(
        output_path(
            f'ld_decay/test/mut_ex/{cell_type}/{chromosome}/{cell_type}_{chromosome}_{gene}_summ_stats.tsv',
            'analysis',
        ),
        sep='\t',
        index=False,
    )


@click.option('--snp-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/tob_freeze_1/bgzip_tabix/v4')
@click.option('--str-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1-chr-specific')
@click.command()
def main(snp_vcf_dir, str_vcf_dir):
    """
    Get all genes that have a lead signal that is not a SNP.
    """
    b = get_batch(name='Lead TR, SNV, and proxy SNV data extraction')

    celltypes = 'gdT,B_intermediate,ILC,Plasmablast'
    #dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive'
    celltypes = celltypes.split(',')
    for cell_type in celltypes:
        for chrom in range(1, 22):
        #for chrom in [22]:
            gene_file = f'gs://cpg-bioheart-test/str/associatr/input_files/240_libraries_tenk10kp1_v2/scRNA_gene_lists/1_min_pct_cells_expressed/{cell_type}/chr{chrom}_{cell_type}_gene_list.json'
            with to_path(gene_file).open() as file:
                genes = json.load(file)
            for gene in genes:

                j = b.new_python_job(
                    name=f'Get pvals/LD of lead variant and closest SNP proxy {cell_type}: {chrom}, {gene}',
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
                    gene
                )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments