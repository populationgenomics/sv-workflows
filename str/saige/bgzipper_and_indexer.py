#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda

"""
This script bgzips and indexes (.csi) VCF files to be used in SAIGE-QTL.
Input MT should be the output of qc_filters_saige.py

 analysis-runner --dataset "bioheart" \
    --description "Hail QC for SAIGE-QTL" \
    --access-level "test" \
    --output-dir "str/saige-qtl/input_files" \
    bgzipper_and_indexer.py --vcf-dir=gs://cpg-bioheart-test/str/saige-qtl/input_files/vcf/v1-chr-specific-rescaled

"""
import logging

import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch


def add_remove_chr_and_index_job(vcf_path):
    """
    Reads a VCF file, it creates an .csi index file and
    and removes "chr" from the chromosome values

    copies and pasted from https://github.com/populationgenomics/saige-tenk10k/blob/main/get_genotype_vcf.py

    Args:
    vcf_path: input path
    """
    # remove chr & add index file using bcftools
    vcf_input = get_batch().read_input(vcf_path)

    vcf_size = to_path(vcf_path).stat().st_size
    storage_required = ((vcf_size // 1024**3) or 1) * 2.2
    bcftools_job = get_batch().new_job(name='remove chr and index vcf')
    bcftools_job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.csi': '{root}.vcf.bgz.csi',
        },
    )
    bcftools_job.image(get_config()['images']['bcftools'])
    bcftools_job.cpu(4)
    bcftools_job.storage(f'{storage_required}Gi')
    # now remove "chr" from chromosome names using bcftools
    bcftools_job.command('for num in {1..22} ; do echo "chr${num} ${num}" >> chr_update.txt; done')
    bcftools_job.command(
        f"""
        bcftools annotate --rename-chrs chr_update.txt {vcf_input} | \\
        bgzip -c > {bcftools_job.output['vcf.bgz']}
        bcftools index -c {bcftools_job.output['vcf.bgz']}
    """,
    )
    logging.info('CV VCF rename/index jobs scheduled!')

    # save both output files
    get_batch().write_output(bcftools_job.output, vcf_path.removesuffix('.vcf.bgz'))


@click.option('--vcf-dir', help='GCS file path to vcf files', type=str)
@click.command()
def main(vcf_dir):
    """
    Runner to bgzip and tabix
    """

    # for chr_index in range(22):  # iterate over chr1-22
    for chr_index in [21]:
        vcf_path = f'{vcf_dir}/hail_filtered_chr{chr_index+1}.vcf.bgz'
        add_remove_chr_and_index_job(vcf_path)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
