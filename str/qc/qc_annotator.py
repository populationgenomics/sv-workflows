#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script annotates the ExpansionHunter MT with the following annotations:
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
- binom_hwep: binomial Hardy-Weinberg equilibrium p-value
- obs_het: proportion of observed heterozygous calls per locus

analysis-runner --access-level "test" --dataset "bioheart" --description "QC annotator" --output-dir "str/polymorphic_run_n2045/annotated_mt/v2" qc_annotator.py \
--mt-path=gs://cpg-bioheart-test/str/polymorphic_run_n2045/mt/v1/str.mt

"""

import hail as hl
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import init_batch, output_path

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

    # sample_qc and variant_qc function
    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)

    # change key (as two variants can have the same locus coordinates, but different REPID)
    mt = mt.annotate_rows(REPID=mt.info.REPID)
    mt = mt.key_rows_by('REPID')

    # clean up repeat length representation - remove 'STR'
    mt = mt.annotate_rows(
        rep_length_alleles=mt.alleles.map(
            lambda x: hl.if_else(
                x.contains('STR'),
                hl.int(x.replace('^<STR([0-9]+)>', '$1')),
                mt.info.REF,
            )
        )
    )
    # motif length
    mt = mt.annotate_rows(motif_length=hl.len(mt.info.RU))

    # bp length array annotation: motif length*repeat length
    mt = mt.annotate_rows(
        bp_length_alleles=mt.rep_length_alleles.map(lambda x: x * mt.motif_length)
    )
    # break up allele_1 and allele_2 into separate columns
    mt = mt.annotate_entries(allele_1_rep_length=hl.int(mt.REPCN.split('/')[0]))
    mt = mt.annotate_entries(
        allele_2_rep_length=hl.if_else(
            hl.len(mt.REPCN.split('/')) == 2,
            hl.int(mt.REPCN.split('/')[1]),
            hl.missing('int32'),
        )
    )
    mt = mt.annotate_entries(
        allele_1_bp_length=mt.allele_1_rep_length * mt.motif_length
    )
    mt = mt.annotate_entries(
        allele_2_bp_length=mt.allele_2_rep_length * mt.motif_length
    )

    # annotate rows with counts of alleles that appear in each VARID
    ht = mt.select_rows(
        alleles_rep_lengths=hl.agg.collect(mt.allele_1_rep_length).extend(
            hl.agg.collect(mt.allele_2_rep_length)
        )
    ).rows()

    # explode the allele_array to create one row for each element in the array
    exploded_table = ht.explode(ht.alleles_rep_lengths)

    # aggregate the counts for each distinct value in the exploded allele_array
    aggregated_table = exploded_table.group_by(exploded_table.REPID).aggregate(
        allele_array_counts=hl.agg.counter(exploded_table.alleles_rep_lengths)
    )

    # finds the mode allele at each locus
    max_key_expr = hl.bind(
        lambda d: hl.bind(
            lambda kv: kv[0], hl.sorted(d.items(), key=lambda kv: -kv[1])[0]
        ),
        aggregated_table.allele_array_counts,
    )
    aggregated_table = aggregated_table.annotate(mode_allele=max_key_expr)

    # annotate the mode allele and dictionary of frequency of each allele into mt
    mt = mt.annotate_rows(aggregated_info=aggregated_table[mt.REPID])
    # number of distinct alleles identified at each locus
    mt = mt.annotate_rows(
        num_alleles=hl.len(mt.aggregated_info.allele_array_counts.values())
    )

    # the difference between each allele and the mode allele
    mt = mt.annotate_entries(
        allele_1_minus_mode=mt.allele_1_rep_length - mt.aggregated_info.mode_allele
    )
    mt = mt.annotate_entries(
        allele_2_minus_mode=mt.allele_2_rep_length - mt.aggregated_info.mode_allele
    )

    # annotate the sum of alleles that are not the mode allele
    mt = mt.annotate_entries(
        allele_1_is_non_mode=hl.if_else(mt.allele_1_minus_mode != 0, True, False)
    )
    mt = mt.annotate_entries(
        allele_2_is_non_mode=hl.if_else(mt.allele_2_minus_mode != 0, True, False)
    )
    mt = mt.annotate_rows(
        sum_allele_1_is_not_mode=hl.agg.sum(hl.cond(mt.allele_1_is_non_mode, 1, 0)),
        sum_allele_2_is_not_mode=hl.agg.sum(hl.cond(mt.allele_2_is_non_mode, 1, 0)),
    )
    mt = mt.annotate_rows(
        sum_alleles_is_not_mode=mt.sum_allele_1_is_not_mode
        + mt.sum_allele_2_is_not_mode
    )
    mt = mt.drop('sum_allele_1_is_not_mode', 'sum_allele_2_is_not_mode')

    # proportion of alleles that are not the mode allele
    mt = mt.annotate_rows(
        prop_alleles_is_not_mode=mt.sum_alleles_is_not_mode
        / hl.sum(mt.aggregated_info.allele_array_counts.values())
    )

    # add in Gymrek binomial HWEP implementation (Hail Query's HWEP implementation only works for biallelic loci)
    # https://github.com/gymrek-lab/TRTools/blob/master/trtools/utils/utils.py#L325
    mt = mt.annotate_rows(n_hom=mt.variant_qc.n_called - mt.variant_qc.n_het)
    mt = mt.annotate_rows(exp_hom_frac=hl.sum(mt.variant_qc.AF**2))
    mt = mt.annotate_rows(n_hom=hl.int32(mt.n_hom))
    mt = mt.annotate_rows(n_called=hl.int32(mt.variant_qc.n_called))
    mt = mt.annotate_rows(
        binom_hwep=hl.binom_test(mt.n_hom, mt.n_called, mt.exp_hom_frac, 'two-sided')
    )
    mt = mt.drop('n_hom', 'n_called', 'exp_hom_frac')

    # proportion of observed heterozygous calls per locus
    mt = mt.annotate_rows(obs_het=mt.variant_qc.n_het / mt.variant_qc.n_called)

    # average 'locus coverage' per locus
    mt = mt.annotate_rows(variant_lc=hl.agg.sum(mt.LC) / mt.variant_qc.n_called)

    # confidence interval (CI) annotations
    mt = mt.annotate_entries(
        allele_1_REPCI=hl.str(mt.REPCI.split('/')[0]),
        allele_2_REPCI=hl.if_else(
            hl.len(mt.REPCN.split('/')) == 2,
            hl.str(mt.REPCI.split('/')[1]),
            hl.missing('str'),
        ),
    )
    mt = mt.annotate_entries(
        allele_1_REPCI_1=hl.int32(mt.allele_1_REPCI.split('-')[0]),
        allele_1_REPCI_2=hl.int32(mt.allele_1_REPCI.split('-')[1]),
    )
    mt = mt.annotate_entries(
        allele_2_REPCI_1=hl.if_else(
            hl.len(mt.allele_2_REPCI.split('-')) == 2,
            hl.int(mt.allele_2_REPCI.split('-')[0]),
            hl.missing('int32'),
        )
    )
    mt = mt.annotate_entries(
        allele_2_REPCI_2=hl.if_else(
            hl.len(mt.allele_2_REPCI.split('-')) == 2,
            hl.int(mt.allele_2_REPCI.split('-')[1]),
            hl.missing('int32'),
        )
    )
    mt = mt.annotate_entries(
        allele_1_REPCI_size=hl.int32(mt.allele_1_REPCI_2 - mt.allele_1_REPCI_1),
        allele_2_REPCI_size=hl.int32(mt.allele_2_REPCI_2 - mt.allele_2_REPCI_1),
    )
    mt = mt.annotate_entries(
        allele_1_REPCI_over_CN=hl.float64(
            mt.allele_1_REPCI_size / mt.allele_1_rep_length
        ),
        allele_2_REPCI_over_CN=hl.float64(
            mt.allele_2_REPCI_size / mt.allele_2_rep_length
        ),
    )
    mt = mt.drop(
        'allele_1_REPCI_1',
        'allele_1_REPCI_2',
        'allele_2_REPCI_1',
        'allele_2_REPCI_2',
        'allele_1_REPCI',
        'allele_2_REPCI',
    )

    # write out mt to GCS path
    mt.write(output_path('str_annotated.mt'))

    # print mt schema
    mt.describe()


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
