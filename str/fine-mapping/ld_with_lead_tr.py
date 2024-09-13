#!/usr/bin/env python3

"""

This script calculates
pairwise LD (Pearson correlation) between the lead TR  and all other variants in the cis window.


Workflow:
1) Extract the genotypes for the lead TR in the cis window .

analysis-runner --dataset "bioheart" \
    --description "Calculate LD for fine-mapped eSTRs with GWAS variants" \
    --access-level "test" \
    --memory='4G' \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr/ld_with_lead_tr" \
    ld_with_lead_tr.py \
    --cell_type CD4_Naive \
    --gene_ensg ENSG00000156475 \
    --chrom chr5


"""


import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def ld_parser(
    cell_type,
    gene_ensg,
    chrom,
    snp_input,
    str_input,
):
    import numpy as np
    from cyvcf2 import VCF

    meta_results = pd.read_csv(
        f'gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results/{cell_type}/{chrom}/{gene_ensg}_100000bp_meta_results.tsv',
        sep='\t',
    )
    meta_results['variant_type'] = np.where(meta_results['motif'].str.contains('-'), 'SNV', 'TR')
    meta_results_tr = meta_results[meta_results['variant_type'] == 'TR']  # filter for TRs
    lead_tr = meta_results_tr[meta_results_tr['pval_meta'] == meta_results_tr['pval_meta'].min()]  # get the lead tr
    lead_tr_coord = str(lead_tr.iloc[0]['pos'])
    lead_tr_motif = lead_tr.iloc[0]['motif']
    lead_str_locus = f'{chrom}_{lead_tr_coord}_{lead_tr_motif}'
    meta_results_max_pos = meta_results['pos'].max()
    meta_results_min_pos = meta_results['pos'].min()
    cis_window = f'{chrom}:{meta_results_min_pos}-{meta_results_max_pos}'

    # Extract the genotypes for SNVs in the cis window
    snp_df = pd.DataFrame(columns=['individual'])
    snp_vcf = VCF(snp_input['vcf'])
    snp_df['individual'] = snp_vcf.samples

    print('Starting to subset SNP VCF for window...')
    for variant in snp_vcf(cis_window):
        gt = variant.gt_types  # extracts GTs as a numpy array
        gt[gt == 3] = 2
        motif = str(variant.INFO.get('RU')).replace('-', '_')
        snp = variant.CHROM + '_' + str(variant.POS) + '_' + motif
        df_to_append = pd.DataFrame(gt, columns=[snp])  # creates a temp df to store the GTs for one locus
        snp_df = pd.concat([snp_df, df_to_append], axis=1)
    print('Finished reading SNP VCF')

    # Extract the genotypes for STRs in the cis window
    str_df = pd.DataFrame(columns=['individual'])
    str_vcf = VCF(str_input['vcf'])
    str_df['individual'] = str_vcf.samples
    print('Starting to subset STR VCF for window...')
    for variant in str_vcf(cis_window):
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

    # merge the STR and SNP GT dfs together
    merged_df = str_df.merge(snp_df, on='individual')

    # calculate pairwise LD with the lead STR
    correlation_series = merged_df.drop(columns='individual').corrwith(merged_df[lead_str_locus])
    correlation_df = pd.DataFrame(correlation_series, columns=['correlation'])
    correlation_df['locus'] = correlation_df.index

    # write results as a tsv file to gcp
    correlation_df.to_csv(
        to_path(
            f'gs://cpg-bioheart-test-analysis/str/associatr/ld_with_lead_tr/{cell_type}/{chrom}/{gene_ensg}_corr.tsv',
        ),
        sep='\t',
        index=False,
    )


@click.option('--snp-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/tob_freeze_1/bgzip_tabix/v4')
@click.option('--str-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1-chr-specific')
@click.option('--cell-type', required=True, help='Cell type')
@click.option('--gene-ensg', required=True, help='Gene ENSG')
@click.option('--chrom', required=True, help='Chromosome')
@click.command()
def main(cell_type, gene_ensg, chrom, snp_vcf_dir, str_vcf_dir):
    b = get_batch(name='Calculate LD for all variants with the lead TR')

    ld_job = b.new_python_job(
        f'LD calc for {chrom}; {cell_type} - {gene_ensg}',
    )

    snp_vcf_path = f'{snp_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'
    str_vcf_path = f'{str_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'

    snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_path, 'tbi': snp_vcf_path + '.tbi'})
    str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'})

    ld_job.call(
        ld_parser,
        cell_type,
        gene_ensg,
        chrom,
        snp_input,
        str_input,
    )

    b.run(wait=False)


if __name__ == '__main__':
    main()
