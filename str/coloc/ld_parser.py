#!/usr/bin/env python3

"""
This script calculates pairwise LD (correlation coefficient) between one STR locus and all SNPs in a specified window.

Outputs the locus-level LD results as a TSV file.

analysis-runner --dataset "bioheart" \
    --description "Calculate LD between STR and SNPs" \
    --access-level "test" \
    --output-dir "str/ld" \
    ld_parser.py --snp-vcf-path=gs://cpg-bioheart-test/saige-qtl/input_files/genotypes/vds1-0/chr22_common_variants.vcf.bgz \
    --str-vcf-path=gs://cpg-bioheart-test/saige-qtl/input_files/genotypes/vds1-0/chr22_common_variants.vcf.bgz \
    --str-locus=22:10510354 \
    --window=22:10510212-10514502 \
    --output-file=ld_results.csv

"""

import click

from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch


def ld_parser(snp_vcf_path: str, str_vcf_path: str, str_locus: str, window: str, output_path: str):

    import pandas as pd
    from cyvcf2 import VCF
    from cpg_utils import to_path

    # copy SNP VCF local because cyVCF2 can only read from a local file
    local_file = 'local.vcf.bgz'
    gcp_file = to_path(snp_vcf_path)
    gcp_file_index = to_path(snp_vcf_path + '.csi')
    gcp_file.copy(local_file)
    gcp_file_index.copy(local_file + '.csi')

    # create empty DF to store the relevant GTs (SNPs)
    df = pd.DataFrame(columns=['individual'])

    # cyVCF2 reads the SNP VCF
    vcf = VCF(local_file)
    for variant in vcf(window):
        geno = variant.gt_types  # extracts GTs as a numpy array
        locus = variant.CHROM + ':' + str(variant.POS)
        df_to_append = pd.DataFrame(geno, columns=[locus])  # creates a temp df to store the GTs for one locus

        # concatenate results to the main df
        df = pd.concat([df, df_to_append], axis=1)

    # extract GTs for the one STR
    # start by copying STR VCF to local
    local_str_file = 'local_str.vcf.bgz'
    gcp_str_file = to_path(str_vcf_path)
    gcp_str_file_index = to_path(str_vcf_path + '.csi')
    gcp_str_file.copy(local_str_file)
    gcp_str_file_index.copy(local_str_file + '.csi')

    # cyVCF2 reads the STR VCF
    str_vcf = VCF(local_str_file)
    for variant in str_vcf(str_locus):
        str_geno = variant.gt_types
        target_data = {'individual': str_vcf.samples, str_locus: str_geno}
        target_df = pd.DataFrame(target_data)

    # merge the two dataframes
    merged_df = df.merge(target_df, on='individual')

    # calculate pairwise correlation of every SNP locus with target STR locus
    correlation_series = merged_df.drop(columns='individual').corrwith(merged_df[str_locus])

    correlation_df = pd.DataFrame(correlation_series, columns=[f'{str_locus}_correlation'])
    correlation_df['locus'] = correlation_df.index
    correlation_df.to_csv(output_path, index=False)


@click.option(
    '--snp-vcf-path',
    help='GCS file path to SNP VCF file.',
    type=str,
)
@click.option(
    '--str-vcf-path',
    help='GCS file path to STR VCF file.',
    type=str,
)
@click.option(
    '--str-locus',
    help='STR locus to calculate LD with.',
    type=str,
)
@click.option(
    '--window',
    help='Window size to calculate LD within.',
    type=str,
)
@click.option(
    '--output-file',
    help='GCS file path to output CSV file.',
    type=str,
)
@click.command()
def main(snp_vcf_path: str, str_vcf_path: str, str_locus: str, window: str, output_file: str):
    b = get_batch('Calculate LD between STR and SNPs')
    ld_job = b.new_python_job(name='LD calculation')
    ld_job.call(ld_parser, snp_vcf_path, str_vcf_path, str_locus, window, output_path(output_file))
    b.run(wait=False)


if __name__ == '__main__':
    main()
