#!/usr/bin/env python3

"""
This script subsets a VCF for samples specified in a text file.

analysis-runner --access-level test --dataset bioheart --description  \
    'VCF sample subsetter'  --output-dir 'str/associatr/tester/vcf' \
    vcf_sample_subsetter.py \
    --input-dir=gs://cpg-bioheart-test/str/dummy_snp_vcf \
    --chromosomes=20 \
    --sample-id-file=gs://cpg-bioheart-test/str/dummy_snp_vcf/tob_ids_str_test.txt
"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


@click.command()
@click.option('--input-dir', help='Input directory for VCF files')
@click.option(
    '--chromosomes',
    default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22',
    help='Chromosome number',
)
@click.option('--sample-id-file', help='File containing CPG IDs to subset')
@click.option('--job-memory', default='4G', help='Job memory')
@click.option('--job-storage', default='10G', help='Job storage')
@click.option('--job-cpu', default=1, help='Job CPU')
def main(input_dir, chromosomes, job_memory, job_storage, job_cpu, sample_id_file):
    """
    Subset a VCF file for samples specified in a text file

    """
    b = get_batch()
    for chrom in chromosomes.split(','):
        input_file = f'{input_dir}/chr{chrom}_common_variants.vcf.bgz'
        vcf_input = b.read_input(input_file)
        sample_id_file = b.read_input(sample_id_file)
        input_file_name = (input_file.split('/')[-1]).split('.')[0]

        bcftools_job = b.new_job(name=f'{input_file_name} Subsetting for Specified Samples')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.memory(job_memory)
        bcftools_job.storage(job_storage)
        bcftools_job.cpu(job_cpu)

        bcftools_job.declare_resource_group(
            vcf_output={
                'vcf.bgz': '{root}.vcf.bgz',
            },
        )

        bcftools_job.command(
            f"""
            bcftools view -S {sample_id_file} -o {bcftools_job.vcf_output['vcf.bgz']} {vcf_input}

            """,
        )

        b.write_output(bcftools_job.vcf_output, (output_path(f'{input_file_name}')))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
