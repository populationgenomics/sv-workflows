#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This Hail Query script subsets a Hail Matrix Table (selects 5 individuals/columns, all rows).

 analysis-runner --dataset "bioheart" \
    --description "subset mt" \
    --access-level "test" \
    --output-dir "str" \
    mt_subset.py --file-path=gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt

"""

import hail as hl
import click

from cpg_utils.hail_batch import output_path, init_batch


@click.option(
    '--file-path',
    help='GCS file path to Hail STR matrix Table',
    type=str,
)
@click.command()
def main(file_path):
    """
    Performs PCA of STR summed repeat length genotypes

    """
    init_batch(worker_memory='highmem')
    mt = hl.read_matrix_table(file_path)

    mt = mt.head(n_rows=15, n_cols=5)
    mt.write(output_path('polymorphic_run_n5_15variants/str_subset.mt'), overwrite=True)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
