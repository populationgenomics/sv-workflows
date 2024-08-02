#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda
"""

This script selects for the 25 samples with LR PacBio data and pulls out the rep lengths.

analysis-runner --access-level test --description "Extract rep lengths for 25 samples with LR PacBio data" --dataset bioheart \
--output-dir potato mt_maker_for_trgt.py


"""

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def main():
    init_batch(worker_memory='highmem')

    mt = hl.read_matrix_table(
        'gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt',
    )

    # Pull out the 25 samples with LR PacBio data
    samples = hl.import_table(
        'gs://cpg-bioheart-test/str/wgs_genotyping/trgt/mapper_trgt.csv',
        delimiter=',',
        key='FS_ID',
        quote=None,
    )
    cpg_ids = samples.s.collect()
    mt = mt.filter_cols(hl.literal(cpg_ids).contains(mt.s))

    # Pull out the rep lengths stored as allele_1_rep_length and allele_2_rep_length
    #mt = mt.annotate_entries(sr_summed_rep_length=mt.allele_1_rep_length + mt.allele_2_rep_length)
    # extract as a table
    #mt.select_entries(mt.sr_summed_rep_length).write(
    #    'gs://cpg-bioheart-test/str/wgs_genotyping/trgt/sr_rep_lengths.mt', overwrite=True,
    #)
    mt = mt.annotate_rows(mode_allele = mt.aggregated_info.mode_allele)
    mt.rows().export(
        'gs://cpg-bioheart-test/str/wgs_genotyping/trgt/analysis-work/mode_allele.tsv.bgz',
    )


if __name__ == '__main__':
    main()
