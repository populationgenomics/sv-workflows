#!/usr/bin/env python3
"""
This Hail Query script outputs a Hail matrix table  from a BGZIPPED VCF file.
 analysis-runner --dataset "bioheart" \
    --description "mt_extractor" \
    --access-level "test" \
    --output-dir "str/polymorphic_run/mt/v1" \
    mt_extractor.py --file-path=gs://cpg-bioheart-test/str/polymorphic_run/bzip/v1/tester_file.vcf.bgz

"""

import click

import hail as hl

from cpg_utils.hail_batch import init_batch, output_path


@click.option(
    '--file-path',
    help='GCS file path to BGZIP VCF file.',
)
@click.command()
def main(file_path):
    """writes a BGZIP VCF as a Hail Matrix Table to a GCS bucket"""

    init_batch(worker_memory='highmem')
    gcs_path = output_path('str.mt')
    hl.import_vcf(file_path, force_bgz=True).write(gcs_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
