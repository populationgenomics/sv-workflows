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
    --output-dir "str/associatr/final-freeze/input_files" \
    qc_filters_associatr.py --mt-path=gs://cpg-bioheart-test/str/associatr/final-freeze/input_files/mt/hail_filtered.mt \
    --version=v1-chr-specific

"""

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch, output_path
config = get_config()




@click.option(
    '--mt-path',
    help='GCS file path to mt (output of qc_annotator.py)',
    type=str,
)

@click.option('--version', help='version of the output files', type=str, default='v1')
@click.command()
def main(
    mt_path,
    version,
):
    """
    Runner to apply QC filters to input MT, and bgzip and tabix.
    """
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



    # wrangle mt, prepare for export_vcf():

    # drop unneccessary columns prior to writing out
    mt = mt.drop(
        'allele_1_rep_length',
        'allele_2_rep_length',
        'allele_1_bp_length',
        'allele_2_bp_length',
        'allele_1_REPCI_over_CN',
        'allele_2_REPCI_over_CN',
        'allele_1_minus_mode',
        'allele_2_minus_mode',
        'allele_1_REPCI_size',
        'allele_2_REPCI_size',
        'allele_1_is_non_mode',
        'allele_2_is_non_mode',
    )

    # Set the QUAL field to missing for all rows
    mt = mt.annotate_entries(QUAL=hl.missing('float64'))

    # Set 'rsid' to REPID
    mt = mt.annotate_rows(rsid=mt.REPID)

    # Key by locus and alleles
    mt = mt.key_rows_by('locus', 'alleles')

    # Remove 'CPG' prefix from ID (associaTR expects ID to be numeric)
    mt = mt.annotate_cols(s_no_prefix=mt.s[3:])
    mt = mt.key_cols_by(s=mt.s_no_prefix)
    mt.drop('s_no_prefix')

    mt.write(
        output_path(f'mt/hail_filtered.mt'),
        overwrite=True,
    )

    for chr_index in range(22):  # iterate over chr1-22
        mt_chr = mt.filter_rows(mt.locus.contig == f'chr{chr_index + 1}')
        gcs_output_path = output_path(f'vcf/{version}/hail_filtered_chr{chr_index+1}.vcf.bgz')
        # needs STR VCF header text to be recognised by associaTR as an ExpansionHunter VCF
        hl.export_vcf(
            mt_chr,
            gcs_output_path,
            append_to_header='gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/hail/STR_header.txt',
            tabix=True,
        )


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
