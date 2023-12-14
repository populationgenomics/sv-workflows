#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script pulls out the locus coordinates of every indel recorded in an existing Hail Matrix Table.

 analysis-runner --dataset "tob-wgs" \
    --description "indel coordinates extractor" \
    --access-level "test" \
    --output-dir "hoptan-str/indel-experiment/" \
    indel_locus_extractor.py --file-path=gs://cpg-tob-wgs-test/mt/v7.mt

"""

import hail as hl
import click

from cpg_utils.hail_batch import output_path, init_batch


def indel_coordinate_extractor(file_path, gcs_path):
    init_batch()
    mt = hl.read_matrix_table(file_path)
    # Filter out monomorphic loci
    mt_filtered = mt.filter_rows(hl.len(mt.alleles) > 1)

    # Filter for indels only
    mt_filtered = mt_filtered.filter_rows(
        hl.is_indel(mt_filtered.alleles[0], mt_filtered.alleles[1])
    )

    # output writing
    mt_filtered.locus.export(gcs_path)


@click.option(
    '--file-path',
    help='GCS file path to Hail matrix table with indels and SNPs',
    type=str,
)
@click.command()
def main(file_path):
    gcs_output_path = output_path(f'indels.tsv', 'analysis')
    indel_coordinate_extractor(file_path, gcs_output_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
