#!/usr/bin/env python3
"""
This Hail Query script prints the dimensions of a VDS
 analysis-runner --dataset "bioheart" \
    --description "mt_extractor" \
    --access-level "test" \
    --output-dir "str/polymorphic_run/mt/v1" \
    vds_counter.py --file-path=gs://cpg-bioheart-test/vds/bioheart1-0.vds

"""

import click

import hail as hl

from cpg_utils.hail_batch import init_batch


@click.option(
    '--file-path',
    help='GCS file path to VDS.',
)
@click.command()
def main(file_path):
    """Print MT dimensons of a VDS"""

    init_batch(worker_memory='highmem')
    vds = hl.vds.read_vds(file_path)
    reference_data_mt = vds.reference_data
    print(f' VDS reference data dimensions: {reference_data_mt.count()}')

    variant_data_mt = vds.variant_data
    print(f' VDS variant data dimensions: {variant_data_mt.count()}')


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
