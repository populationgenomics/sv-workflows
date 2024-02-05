#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This Hail Query script joins two Hail matrix tables together to form one matrix table, outputting to a GCS bucket.
The rows of the matrix table passed in as mt-path-1 will be the rows for the final mt.
The columns of both matrix tables will be in the final mt.

 analysis-runner --dataset "bioheart" \
    --description "mt_joiner" \
    --access-level "test" \
    --output-dir "str/polymorphic_run/combined_mt/v1" \
    mt_joiner.py --mt-path-1=gs://cpg-bioheart-test/str/polymorphic_run/bzip/v1/tester_file.mt \
    --mt-path-2=gs://cpg-bioheart-test/str/polymorphic_run/bzip/v1/tester_file.mt
"""

import hail as hl
import click

from cpg_utils.hail_batch import output_path, init_batch


@click.option(
    '--mt-path-1',
    help='GCS file path to first Matrix Table',
)
@click.option(
    '--mt-path-2',
    help='GCS file path to second Matrix Table',
)
@click.command()
def main(mt_path_1, mt_path_2):
    """joins two Hail matrix tables together and writes combined mt to GCS"""

    init_batch()
    mt_1 = hl.read_matrix_table(mt_path_1)

    #rekey by REPID
    mt_1 = mt_1.annotate_rows(REPID = mt_1.info.REPID)
    mt_1 = mt_1.key_rows_by(REPID = mt_1['REPID'])

    print(f'{mt_path_1} dimensions: {mt_1.count()}')

    mt_2 = hl.read_matrix_table(mt_path_2)

    #rekey by REPID
    mt_2 = mt_2.annotate_rows(REPID = mt_2.info.REPID)
    mt_2 = mt_2.key_rows_by(REPID = mt_2['REPID'])

    print(f'{mt_path_2} dimensions: {mt_2.count()}')

    mt_joined = mt_1.union_cols(mt_2)
    print(f'Joined matrix table dimensions: {mt_joined.count()}')

    gcs_output_path = output_path(f'combined.mt', 'analysis')
    mt_joined.write(gcs_output_path, overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
