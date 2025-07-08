#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda
"""
This script applies filters to a ExpansionHunter MT, and outputs chromosome-specific VCF ready for input into associaTR.
Input MT should be the output of qc_annotator.py

Applied filters:
- remove monomorphic variants
- locus-level callrate >=0.9
- sample-level call rate >=0.99
- observed heterozygosity >=0.00995 (corresponds to a MAF = 0.5%)
- locus-level HWEP (binom definition) >=10^-6
- set calls outside of [-30,20] relative to mode to NA
- enforce locus-level call rate >=0.9


 analysis-runner --dataset "bioheart" \
    --description "Hail QC for associaTR" \
    --access-level "test" \
    --output-dir "str/polymorphic_run/mt/bioheart_tob/v1_n1925/default_filters" \
    qc_filters_associatr.py --mt-path=gs://cpg-bioheart-test/str/polymorphic_run/mt/bioheart_tob/v1_n1925/str_annotated.mt

"""

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


def qc_filter(mt_path, version):
    """
    Applies QC filters to the input MT
    """
    from cpg_utils.hail_batch import init_batch, output_path

    init_batch(worker_memory='highmem')

    # read in mt
    mt = hl.read_matrix_table(mt_path)

    # remove monomorphic variants, set locus level call rate >=0.9, observed heterozygosity >=0.00995, locus level HWEP (binom definition) >=10^-6
    mt = mt.filter_rows(
        (mt.num_alleles > 1) & (mt.variant_qc.call_rate >= 0.9) & (mt.obs_het >= 0.00995) & (mt.binom_hwep >= 0.000001),
    )

    # set sample level call rate >=0.99
    mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.99)

    # set big expansions/deletions beyond [-30,20] relative to mode allele to NA
    condition_allele_1 = (mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)
    mt = mt.annotate_entries(
        GT=hl.if_else(
            condition_allele_1,
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            condition_allele_1,
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.SO,
        ),
    )
    condition_allele_2 = (mt.allele_2_minus_mode > 20) | (mt.allele_2_minus_mode < -30)
    mt = mt.annotate_entries(
        GT=hl.if_else(
            condition_allele_2,
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            condition_allele_2,
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.SO,
        ),
    )

    # calculate proportion of GTs that are defined per locus (after applying call-level filters, variant_qc.call_rate is not accurate anymore)
    mt = mt.annotate_rows(
        prop_GT_exists=hl.agg.count_where(hl.is_defined(mt.GT)) / (mt.variant_qc.n_called + mt.variant_qc.n_not_called),
    )
    # re-enforce locus level call rate >=0.9
    mt = mt.filter_rows(mt.prop_GT_exists >= 0.9)

    # sum the alleles
    mt = mt.annotate_entries(
    sum_alleles = mt.allele_1_rep_length + mt.allele_2_rep_length
    )
    # get the minimum dosage per locus
    mt = mt.annotate_rows(
        min_dosage = hl.agg.filter(
            hl.is_defined(mt.sum_alleles),
            hl.agg.min(mt.sum_alleles)
        )
    )



    # Set the QUAL field to missing for all rows
    mt = mt.annotate_entries(QUAL=hl.missing('float64'))

    # Set 'rsid' to REPID
    mt = mt.annotate_rows(rsid=mt.REPID)

    # Key by locus and alleles
    mt = mt.key_rows_by('locus', 'alleles')

    mt.write(output_path('str_annotated.mt'))

    mt.rows().export(output_path('str_annotated_rows.tsv.bgz'))

@click.option(
    '--mt-path',
    help='GCS file path to mt (output of qc_annotator.py)',
    type=str,
)
@click.option('--hail-storage', help='Hail storage', type=str, default='0G')
@click.option('--hail-cpu', help='Hail CPU', type=int, default=1)
@click.option('--hail-memory', help='Hail memory', type=str, default='standard')
@click.option('--version', help='version of the output files', type=str, default='v1')
@click.command()
def main(
    mt_path,
    hail_storage,
    hail_cpu,
    hail_memory,
    version,
):
    """
    Runner to apply QC filters to input MT, and bgzip and tabix.
    """

    b = get_batch('Apply QC filters to MT and export chr-specific VCFs')

    hail_job = b.new_python_job(name='QC filters')
    hail_job.image(config['workflow']['driver_image'])
    hail_job.storage(hail_storage)
    hail_job.cpu(hail_cpu)
    hail_job.memory(hail_memory)
    hail_job.call(qc_filter, mt_path, version)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
