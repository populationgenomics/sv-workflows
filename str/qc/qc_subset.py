#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""


analysis-runner --access-level "test" --dataset "bioheart" --description "QC annotator" --output-dir "str/polymorphic_run/mt/bioheart_tob/v1_n2412/v1-default-filters/v1_n1925" qc_subset.py \
--mt-path=gs://cpg-bioheart-test/str/polymorphic_run/mt/bioheart_tob/v1_n2412/v1-default-filters/str_filtered.mt

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
    print(f'MT dimensions: {mt.count()}')

    bioheart_ids = pd.read_csv('gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/bioheart_n975_sample_covariates.csv')['sample_id']
    tob_ids = pd.read_csv('gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/tob_n950/covariates/6_rna_pcs/CD4_TCM_covariates.csv')['sample_id']
    samples = tob_ids.to_list() + bioheart_ids.to_list()

    # filter the MT to only include samples in the sample list
    mt = mt.filter_cols(hl.literal(samples).contains(mt.s))


    # write out mt to GCS path
    mt.write(output_path('str_filtered.mt'))

    # print mt schema
    mt.describe()


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
