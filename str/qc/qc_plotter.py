#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script generates plots from the Hail Matrix Table

"""

import hail as hl
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import init_batch

config = get_config()

def main():
    init_batch(worker_memory='highmem')
    mt = hl.read_matrix_table('gs://cpg-bioheart-test/str/polymorphic_run/mt_joiner/n_2045.mt')
    print(f' MT dimensions: {mt.count()}')




if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
