#!/usr/bin/env python3

"""
This script receives a VCF file and performs BGZIP and TABIX

analysis-runner --access-level test --dataset bioheart --description  \
    'VCF combiner'  --output-dir 'tenk10k/str/associatr/final_freeze/methyl_eqtl/input_files' \
    bgzip_tabix.py \
    --input-vcf=gs://gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/methyl_eqtl/input_files/finemapped_etrs.vcf

"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


@click.command()
@click.option('--input-vcf', help='Input directory for VCF files')

@click.option('--job-memory', default='4G', help='Job memory')
@click.option('--job-storage', default='10G', help='Job storage')
@click.option('--job-cpu', default=1, help='Job CPU')
def main(input_vcf, chromosomes, job_memory, job_storage, job_cpu):
    """
    BGZIPs and TABIX a GZIPPED/unzipped VCF input file and writes it to a GCS bucket

    """
    b = get_batch()


    vcf_input = b.read_input(input_vcf)
    input_file_name = (input_vcf.split('/')[-1]).split('.')[0]

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
        tabix -f -p vcf {bcftools_job.vcf_output['vcf.bgz']}  > {bcftools_job.vcf_output['vcf.bgz.tbi']}

        """,
    )
    b.write_output(bcftools_job.vcf_output, (output_path(f'{input_file_name}')))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
