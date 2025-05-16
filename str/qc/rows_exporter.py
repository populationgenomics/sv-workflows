#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script exports the rows of the MT to a TSV file.

analysis-runner --access-level "test" --dataset "bioheart" --description "QC annotator" --output-dir "str/polymorphic_run/mt/bioheart_tob/v1_n1925" rows_exporter.py \
--mt-path=gs://cpg-bioheart-test/str/polymorphic_run/mt/bioheart_tob/v1_n1925/str_annotated.mt

"""

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch, output_path
import pandas as pd

config = get_config()


@click.option('--mt-path', help='GCS Path to the input MT')
@click.command()
def main(mt_path):
    """
    Annotates the ExpansionHunter MT, and outputs annotated MT to GCS
    """

    init_batch(worker_memory='highmem')
    mt = hl.read_matrix_table(mt_path)

    mt.rows().export('gs://cpg-bioheart-test/str/polymorphic_run/mt/bioheart_tob/v1_n1925/rows.tsv.bgz')

    # print mt schema
    mt.describe()


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
