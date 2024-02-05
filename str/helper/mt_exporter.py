#!/usr/bin/env python3

"""
This script receives a gzipped VCF file and creates a Hail Matrix table, written to a GCS bucket.

analysis-runner --access-level test --dataset bioheart --description  \
    'VCF combiner' --memory 32G --storage 50G --output-dir 'str/polymorphic_run/mt/bioheart/v1' \
    mt_exporter.py \
    --input-file=gs://cpg-bioheart-test/str/polymorphic_run/merge_str/tester_file.gz
"""
import hail as hl
import click


from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path, init_batch


config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


def mt_writer(file_path, gcs_path):
    """ writes a BGZIP VCF as a Hail Matrix Table to a GCS bucket """
    init_batch()
    hl.import_vcf(file_path, force_bgz=True).write(gcs_path, overwrite=True)


@click.command()
@click.option('--input-file', help='Parent input file path for gzipped VCF')
@click.option('--job-memory', default='8G', help='Job memory')
@click.option('--job-storage', default='8G', help='Job storage')
def main(input_file, job_memory, job_storage):
    """
    BGZIPs a GZIPPED VCF input file and writes it to a GCS bucket as a Hail Matrix Table.

    """
    b = get_batch()
    vcf_input = b.read_input(input_file)
    bcftools_job = b.new_job(name=f'{input_file} Files prep')
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.memory(job_memory)
    bcftools_job.storage(job_storage)

    bcftools_job.declare_resource_group(
        vcf_output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
        }
    )

    bcftools_job.command(
        f"""

        bcftools view {vcf_input} | bgzip -c >  ${BATCH_TMPDIR}/vcf.bgz

        """
    )

    mt_output_path = to_path(output_path('str.mt', 'analysis'))
    mt_writer_job = b.new_python_job(name='mt_writer job')
    mt_writer_job.call(mt_writer, f'${BATCH_TMPDIR}/vcf.bgz', mt_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
