#!/usr/bin/env python3

"""
This script receives a gzipped VCF file and performs BGZIP, writing files to output (prep files for Hail Query).

analysis-runner --access-level test --dataset bioheart --description  \
    'VCF combiner' --memory 32G --storage 50G --output-dir 'str/polymorphic_run/mt/bioheart/v1' \
    bzip.py \
    --input-file=gs://cpg-bioheart-test/str/polymorphic_run/merge_str/tester_file.gz
"""

import click


from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path


config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


@click.command()
@click.option('--input-file', help='Parent input file path for gzipped VCF')
@click.option('--job-memory', default='8G', help='Job memory')
@click.option('--job-storage', default='8G', help='Job storage')
@click.option('--job-cpu', default=4, help='Job CPU')
def main(input_file, job_memory, job_storage, job_cpu):
    """
    BGZIPs a GZIPPED VCF input file and writes it to a GCS bucket as a Hail Matrix Table.

    """
    b = get_batch()
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
        }
    )

    bcftools_job.command(
        f"""

        bcftools view {vcf_input} | bgzip -c >  {bcftools_job.vcf_output['vcf.bgz']}


        """
    )
    b.write_output(
        bcftools_job.vcf_output, (output_path(f'{input_file_name}', 'analysis'))
    )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
