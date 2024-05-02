#!/usr/bin/env python3
"""
This Hail Query script checks subsetting consistency between multi-cohort VDS and single-cohort VDS.
 analysis-runner --dataset "bioheart" \
    --description "mt_extractor" \
    --access-level "full" \
    --output-dir "str/" \
    vds_checker.py --tenk10k-filepath=gs://cpg-bioheart-main/vds/tenk10k1-0.vds --bioheart-filepath=gs://cpg-bioheart-main/vds/bioheart1-0.vds

"""

import click

import hail as hl

from cpg_utils.hail_batch import init_batch


@click.option(
    '--tenk10k-filepath',
    help='GCS file path to VDS.',
)
@click.option(
    '--bioheart-filepath',
    help='GCS file path to VDS.',
)
@click.option(
    '--num-samples',
    help='Number of samples randomly selected to compare ',
    default=15
)
@click.command()
def main(tenk10k_filepath, bioheart_filepath, num_samples):
    """Check if subsetting VDS is the same as producing VDS from scratch using only subsetted samples"""

    init_batch(worker_memory='highmem')

    tenk10k_vds = hl.vds.read_vds(tenk10k_filepath)
    bioheart_vds = hl.vds.read_vds(bioheart_filepath)

    tenk10k_reference_data_mt = tenk10k_vds.reference_data
    bioheart_reference_data_mt = bioheart_vds.reference_data
    print(f' VDS reference data dimensions: {tenk10k_reference_data_mt.count()}')
    print(f' VDS reference data dimensions: {bioheart_reference_data_mt.count()}')

    tenk10k_variant_data_mt = tenk10k_vds.variant_data
    bioheart_variant_data_mt = bioheart_vds.variant_data
    print(f' VDS variant data dimensions: {tenk10k_variant_data_mt.count()}')
    print(f' VDS variant data dimensions: {bioheart_variant_data_mt.count()}')

    bioheart_only_sgids = bioheart_vds.variant_data.s.collect()
    # randomly sample 15 IDs
    bioheart_sampled_ids = random.sample(bioheart_only_sgids, k=num_samples)
    # filter to just test subset
    filtered_tenk10k_vds = hl.vds.filter_samples(tenk10k_vds, bioheart_sampled_ids)
    filtered_bioheart_vds = hl.vds.filter_samples(bioheart_vds, bioheart_sampled_ids)

    # densify both
    filtered_tenk10k_dense = hl.vds.to_dense_mt(filtered_tenk10k_vds)
    bioheart_dense = hl.vds.to_dense_mt(bioheart_vds)
    print(f' Filtered tenk10K dense_mt data dimensions: {filtered_tenk10k_dense.count()}')
    print(f' BioHEART dense_mt data dimensions: {bioheart_dense.count()}')

    # compare vds represenations
    print(f'Variant data VDS same? {filtered_tenk10k_vds.variant_data._same(bioheart_vds.variant_data)}')

    # compare dense_mt's
    print(f'Dense Matrix Tables same? {filtered_tenk10k_dense._same(bioheart_dense)}')


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
