#!/usr/bin/env python3

"""
This script receives a gzipped VCF file and performs BGZIP, writing files to output (prep files for Hail Query).

analysis-runner --access-level test --dataset bioheart --description  \
    'VCF combiner' --memory 32G --storage 50G --output-dir 'str/associatr/common_variant_snps' \
    bzip.py \
    --input-dir=gs://cpg-bioheart-test/str/associatr/common_variants_snps \
    --chromosomes=20
"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


@click.command()
@click.option('--input-dir', help='Input directory for VCF files')
@click.option(
    '--chromosomes', default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22', help='Chromosome number',
)
@click.option('--job-memory', default='8G', help='Job memory')
@click.option('--job-storage', default='8G', help='Job storage')
@click.option('--job-cpu', default=4, help='Job CPU')
def main(input_dir, chromosomes, job_memory, job_storage, job_cpu):
    """
    BGZIPs and TABIX a GZIPPED/unzipped VCF input file and writes it to a GCS bucket

    """
    b = get_batch()
    for chrom in chromosomes.split(','):
        input_file = f'{input_dir}/hail_filtered_chr{chrom}.vcf'
        vcf_input = b.read_input(input_file)
        input_file_name = (input_file.split('/')[-1]).split('.')[0]

        bcftools_job = b.new_job(name=f'{input_file_name} BGZIPPING')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.memory(job_memory)
        bcftools_job.storage(job_storage)
        bcftools_job.cpu(job_cpu)

        bcftools_job.declare_resource_group(
            vcf_output={
                'vcf.bgz': '{root}.vcf.bgz',
                'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            },
        )

        bcftools_job.command(
            f"""

            bcftools view {vcf_input} | bgzip -c >  {bcftools_job.vcf_output['vcf.bgz']}
            tabix -f -p vcf {bcftools_job.vcf_output['vcf.bgz']}  > {bcftools_job.vcf_output['vcf.gz.tbi']}

            """,
        )
        b.write_output(bcftools_job.vcf_output, (output_path(f'{input_file_name}')))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
