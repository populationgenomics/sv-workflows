#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda
"""
This script applies filters to a ExpansionHunter MT, and outputs chromosome-specific VCF ready for input into SAIGE-QTL.
Input MT should be the output of qc_annotator.py

Applied filters:
- remove monomorphic variants
- locus-level callrate >=0.9
- sample-level call rate >=0.99
- observed heterozygosity >=0.00995 (corresponds to a MAF = 0.5%)
- locus-level HWEP (binom definition) >=10^-6
- set calls outside of [-30,20] relative to mode to NA
- enforce locus-level call rate >=0.9
- set NA alleles to mode allele for the respective locus (SAIGE-QTL does not support missing alleles)


 analysis-runner --dataset "bioheart" \
    --description "Hail QC for SAIGE-QTL" \
    --access-level "test" \
    --output-dir "str/saige-qtl/input_files" \
    qc_filters_saige.py --mt-path=gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt \
    --version=v1-chr-specific

"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch

config = get_config()


def qc_filter(mt_path, version):
    """
    Applies QC filters to the input MT

    """
    from cpg_utils.hail_batch import init_batch,output_path
    import hail as hl


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

    # create DS entry field (summed repeat length across both alleles); set to mode*2 if GT is missing
    mt = mt.annotate_entries(
        DS=hl.if_else(
            hl.is_defined(mt.GT), mt.allele_1_rep_length + mt.allele_2_rep_length, mt.aggregated_info.mode_allele * 2,
        ),
    )

    # wrangle mt, prepare for export_vcf():

    # drop unneccessary entr fields prior to writing out
    mt = mt.drop(
        'GT',
        'ADFL',
        'ADIR',
        'ADSP',
        'LC',
        'REPCI',
        'REPCN',
        'SO',
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
    mt = mt.drop(mt.info)

    # Set the QUAL field to missing for all rows
    mt = mt.annotate_entries(QUAL=hl.missing('float64'))

    # Set 'rsid' to REPID
    mt = mt.annotate_rows(rsid=mt.REPID)

    # Key by locus and alleles (dummy REF = A and ALT = C)
    mt = mt.key_rows_by(locus=mt['locus'], alleles=hl.str(mt['rep_length_alleles']))
    mt = mt.annotate_rows(str_array_field=["A", "C"])
    mt = mt.key_rows_by(locus=mt['locus'], alleles=mt['str_array_field'])

    for chr_index in range(22):  # iterate over chr1-22
        mt_chr = mt.filter_rows(mt.locus.contig == f'chr{chr_index + 1}')
        gcs_output_path = output_path(f'vcf/{version}/hail_filtered_chr{chr_index+1}.vcf.bgz')
        hl.export_vcf(
            mt_chr,
            gcs_output_path,
        )


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

    b = get_batch('Apply QC filters to MT and export chr-specific VCFs')

    hail_job = b.new_python_job(name='QC filters')
    hail_job.image(config['workflow']['driver_image'])
    hail_job.call(qc_filter, mt_path, version)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
