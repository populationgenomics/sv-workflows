#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script outputs a CSV file containing the CPG ID, external ID, and karyotypic sex from the SampleQC Hail Table.

 analysis-runner --dataset "tob-wgs" \
    --description "karyotypic sex extractor" \
    --access-level "test" \
    --output-dir "hoptan-str/sample-sex-mapping/" \
    karyotype_sex_extractor.py --file-path=gs://cpg-bioheart-test/large_cohort/1-0/sample_qc.ht/

"""

import hail as hl
import click

from cpg_utils.hail_batch import output_path, init_batch


def karyotype_sex_extractor(file_path, gcs_path):
    init_batch()
    sample_qc_table = hl.read_table(file_path)

    # select CPG ID, external ID, and sex karyotype (CPG ID is the key, so don't need to explicitly select)
    karyotype_sex_table = sample_qc_table.select(
        sample_qc_table.external_id, sample_qc_table.sex_karyotype
    )

    # output writing
    karyotype_sex_table.export(gcs_path, delimiter=",")


@click.option(
    '--file-path',
    help='GCS file path to SampleQC Hail Table',
    type=str,
)
@click.command()
def main(file_path):
    gcs_output_path = output_path(f'sample_karyotype_sex_mapping.csv', 'analysis')
    karyotype_sex_extractor(file_path, gcs_output_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
