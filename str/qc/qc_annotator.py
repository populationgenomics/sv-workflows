#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script annotates the MT with the following annotations:
- rep_length_alleles: repeat length of each allele
- motif_length: length of the motif
- bp_length_alleles: bp length of each allele
- allele_1_rep_length: repeat length of allele 1
- allele_2_rep_length: repeat length of allele 2
- allele_1_bp_length: bp length of allele 1
- allele_2_bp_length: bp length of allele 2
- aggregated_info.mode_allele:mode allele at each locus
- num_alleles: number of alleles at each locus
- allele_1_minus_mode: difference between allele 1 and mode allele
- allele_2_minus_mode: difference between allele 2 and mode allele
- sum_alleles_is_not_mode: count of non-mode alleles at each locus
- prop_alleles_is_not_mode: proportion of alleles that are not the mode allele per locus

"""

import hail as hl
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import init_batch

config = get_config()

def main():
    init_batch(worker_memory='highmem')
    mt = hl.read_matrix_table('gs://cpg-bioheart-test/str/polymorphic_run/mt_joiner/n_2045.mt')
    print(f'MT dimensions: {mt.count()}')

    ## annotations
    # inspired by Jupyter notebook on 5M-200.mt

    #wrangle repeat length representation - remove 'STR'
    mt = mt.annotate_rows(
        rep_length_alleles=mt.alleles.map(
            lambda x: hl.if_else(
                x.contains("STR"),
                hl.int(x.replace("^<STR([0-9]+)>", "$1")),
                mt.info.REF,
            )
        )
    )
    # motif length
    mt = mt.annotate_rows(motif_length = hl.len(mt.info.RU))

    #bp length annotation: motif length*repeat length
    mt = mt.annotate_rows(
        bp_length_alleles=mt.rep_length_alleles.map(
            lambda x:
                x*mt.motif_length
            )
        )
    #further wrangling of repeat length specs
    mt = mt.annotate_entries(allele_1_rep_length=hl.int(mt.REPCN.split("/")[0]))
    mt = mt.annotate_entries(allele_2_rep_length = hl.if_else(hl.len(mt.REPCN.split("/")) ==2, hl.int(mt.REPCN.split("/")[1]), hl.missing('int32')))
    mt = mt.annotate_entries(allele_1_bp_length= mt.allele_1_rep_length*mt.motif_length)
    mt = mt.annotate_entries(allele_2_bp_length= mt.allele_2_rep_length*mt.motif_length)

    #annotate rows with counts of alleles that appear in each VARID
    ht = mt.select_rows(
        alleles_rep_lengths = hl.agg.collect(mt.allele_1_rep_length)
            .extend(hl.agg.collect(mt.allele_2_rep_length))
    ).rows()

    # Explode the allele_array to create one row for each element in the array
    exploded_table = ht.explode(ht.alleles_rep_lengths)


    # Aggregate the counts for each distinct value in the exploded allele_array
    aggregated_table = exploded_table.group_by(exploded_table.REPID).aggregate(
        allele_array_counts=hl.agg.counter(exploded_table.alleles_rep_lengths)
    )

    #finds the mode allele at each locus
    max_key_expr = hl.bind(lambda d: hl.bind(lambda kv: kv[0], hl.sorted(d.items(), key=lambda kv: -kv[1])[0]), aggregated_table.allele_array_counts)
    aggregated_table = aggregated_table.annotate(mode_allele=max_key_expr)

    mt = mt.annotate_rows(aggregated_info = aggregated_table[mt.REPID])
    mt = mt.annotate_rows(num_alleles = hl.len(mt.aggregated_info.allele_array_counts.values()))

    mt = mt.annotate_entries(allele_1_minus_mode = mt.allele_1_rep_length- mt.aggregated_info.mode_allele)
    mt = mt.annotate_entries(allele_2_minus_mode = mt.allele_2_rep_length-mt.aggregated_info.mode_allele)

    mt = mt.annotate_entries(allele_1_is_non_mode = hl.if_else(mt.allele_1_minus_mode!=0, True, False))
    mt = mt.annotate_entries(allele_2_is_non_mode = hl.if_else(mt.allele_2_minus_mode!=0, True, False))

    mt = mt.annotate_rows(sum_allele_1_is_not_mode=hl.agg.sum(hl.cond(mt.allele_1_is_non_mode, 1, 0)),
                         sum_allele_2_is_not_mode=hl.agg.sum(hl.cond(mt.allele_2_is_non_mode, 1, 0)))
    mt = mt.annotate_rows(sum_alleles_is_not_mode = mt.sum_allele_1_is_not_mode + mt.sum_allele_2_is_not_mode)
    mt = mt.drop("sum_allele_1_is_not_mode")
    mt = mt.drop("sum_allele_2_is_not_mode")
    mt = mt.annotate_rows(prop_alleles_is_not_mode = mt.sum_alleles_is_not_mode/hl.sum(mt.aggregated_info.allele_array_counts.values()))

    #write out mt
    mt.write('gs://cpg-bioheart-test/str/polymorphic_run/mt_joiner/n_2045_annotated.mt')

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
