#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda
"""

This script selects for the 25 samples with LR PacBio data and pulls out the rep lengths.

analysis-runner --access-level test --description "Extract rep lengths for 25 samples with LR PacBio data" --dataset bioheart \
--output-dir potato mt_maker_for_trgt_sr.py


"""

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def main():
    init_batch(worker_memory='highmem')

    mt = hl.import_vcf(
        'gs://cpg-bioheart-test/str/wgs_genotyping/trgt/output_trgt_all_v2.vcf.gz', force=True, array_elements_required=False,
    )
    # Load the mapping file (assuming it's a CSV with columns 'old_id' and 'new_id')
    mapping_table = hl.import_table(
        'gs://cpg-bioheart-test/str/wgs_genotyping/trgt/mapper_trgt.csv',
        delimiter=',',
        key='FS_ID',
        quote=None,
    )

    mt = mt.annotate_cols(new_id=mapping_table[mt.s].s)
    mt = mt.key_cols_by('new_id')
    mt = mt.drop(mt.s)

    # Rename an existing column (e.g., rename 'old_column_name' to 'new_column_name')
    mt = mt.rename({'new_id': 's'})
    mt = mt.annotate_rows(REPID=mt.info.TRID)
    mt = mt.key_rows_by('REPID')
    sr = hl.read_matrix_table('gs://cpg-bioheart-test/str/wgs_genotyping/trgt/sr_rep_lengths.mt')
    mt = mt.annotate_entries(sr_summed_rep_length=sr[mt.row_key, mt.col_key].sr_summed_rep_length)

    mt.checkpoint('gs://cpg-bioheart-test/str/wgs_genotyping/trgt/trgt_sr_25_v2.mt')


    qc_table = hl.import_table(
        'gs://cpg-bioheart-test/str/polymorphic_run/mt/bioheart_tob/v1_n2412/v1-default-filters/str_annotated_rows.tsv.bgz',
        delimiter='\t',
        quote=None,
    )
    qc_table = qc_table.key_by('REPID')
    repid_list = qc_table.aggregate(hl.agg.collect(qc_table.REPID))
    repid_set = hl.set(repid_list)
    mt = mt.filter_rows(repid_set.contains(mt.REPID))

    mt.checkpoint('gs://cpg-bioheart-test/str/wgs_genotyping/trgt/trgt_sr_25_v2_filtered.mt')

    mt = mt.annotate_entries(lr_summed_rep_length=hl.int(mt.MS[0]))
    mt = mt.annotate_entries(strict_concord=hl.if_else(mt.sr_summed_rep_length == mt.lr_summed_rep_length, 1, 0))
    mt = mt.annotate_entries(
        off_by_one_concord=hl.if_else(hl.abs(mt.sr_summed_rep_length - mt.lr_summed_rep_length) <= 1, 1, 0),
    )
    mt = mt.annotate_rows(
        prop_strict_concord=hl.agg.sum(mt.strict_concord) / 25,
        prop_off_by_one_concord=hl.agg.sum(mt.off_by_one_concord) / 25,
    )
    mt.checkpoint('gs://cpg-bioheart-test/str/wgs_genotyping/trgt/trgt_sr_25_v2_filtered_annotated.mt')
    mt.rows().export('gs://cpg-bioheart-test/str/wgs_genotyping/trgt/analysis-work/final-freeze/filtered_trgt_sr_25_v2_rows.tsv.bgz')




if __name__ == '__main__':
    main()
