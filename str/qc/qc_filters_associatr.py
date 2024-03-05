#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda
"""
This script applies filters to a ExpansionHunter MT, and outputs a VCF ready for input into associaTR.
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
    --output-dir "str/associatr/input_files" \
    qc_filters_associatr.py --mt-path=gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt

"""

import hail as hl
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import get_batch, init_batch, output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


def qc_filter(mt_path, gcs_path):
    """
    Applies QC filters to the input MT
    """
    init_batch(worker_memory='highmem')

    # read in mt
    mt = hl.read_matrix_table(mt_path)

    # remove monomorphic variants, set locus level call rate >=0.9, observed heterozygosity >=0.00995, locus level HWEP (binom definition) >=10^-6
    mt = mt.filter_rows(
        (mt.num_alleles > 1)
        & (mt.variant_qc.call_rate >= 0.9)
        & (mt.obs_het >= 0.00995)
        & (mt.binom_hwep >= 0.000001)
    )

    # set sample level call rate >=0.99
    mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.99)

    # set big expansions/deletions beyond [-20,30] relative to mode allele to NA
    mt = mt.annotate_entries(
        GT=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.SO,
        ),
    )

    mt = mt.annotate_entries(
        GT=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            ((mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)),
            hl.missing('str'),
            mt.SO,
        ),
    )

    # calculate proportion of GTs that are defined per locus (after applying call-level filters, variant_qc.call_rate is not accurate anymore)
    mt = mt.annotate_rows(
        prop_GT_exists=hl.agg.count_where(hl.is_defined(mt.GT))
        / (mt.variant_qc.n_called + mt.variant_qc.n_not_called)
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
    # Set 'alleles' to 'rep_length_alleles' (transform to string type to meet export_vcf() specs)
    mt = mt.annotate_rows(alleles=hl.map(lambda x: hl.str(x), mt.rep_length_alleles))

    # Key by locus and alleles
    mt = mt.key_rows_by('locus', 'alleles')

    # Remove 'CPG' prefix from ID (associaTR expects ID to be numeric)
    mt = mt.annotate_cols(s_no_prefix=mt.s[3:])
    mt = mt.key_cols_by(s=mt.s_no_prefix)
    mt.drop('s_no_prefix')

    # needs STR VCF header text to be recognised by associaTR as an ExpansionHunter VCF
    hl.export_vcf(
        mt,
        gcs_path,
        append_to_header='gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/hail/STR_header.txt',
    )


@click.option(
    '--mt-path',
    help='GCS file path to mt (output of qc_annotator.py)',
    type=str,
)
@click.option('--hail-storage', help='Hail storage', type=str, default='20G')
@click.option('--hail-cpu', help='Hail CPU', type=int, default=4)
@click.option('--hail-memory', help='Hail memory', type=str, default='standard')
@click.option('--bcftools-storage', help='BCFTOOLS storage', type=str, default='300G')
@click.option('--bcftools-cpu', help='BCFTOOLS CPU', type=int, default=8)
@click.option('--bcftools-memory', help='BCFTOOLS memory', type=str, default='standard')
@click.option('--version', help='version of the output files', type=str, default='v1')
@click.command()
def main(
    mt_path,
    hail_storage,
    hail_cpu,
    hail_memory,
    bcftools_storage,
    bcftools_cpu,
    bcftools_memory,
    version,
):
    """
    Runner to apply QC filters to input MT, and bgzip and tabix.
    """

    b = get_batch()

    gcs_output_path = output_path(f'vcf/{version}/hail_filtered.vcf')
    #hail_job = b.new_python_job(name=f'QC filters')
    #hail_job.image(config['workflow']['driver_image'])
    #hail_job.storage(hail_storage)
    #hail_job.cpu(hail_cpu)
    #hail_job.memory(hail_memory)
    #hail_job.call(qc_filter, mt_path, gcs_output_path)

    bcftools_job = b.new_job(name=f'bgzip and tabix the Hail output VCF')
    #bcftools_job.depends_on(hail_job)
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.storage(bcftools_storage)
    bcftools_job.cpu(bcftools_cpu)
    bcftools_job.memory(bcftools_memory)
    vcf = b.read_input(gcs_output_path)

    bcftools_job.declare_resource_group(
        vcf_output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    bcftools_job.command(
        f"""
    set -ex;
    echo "Compressing";
    bcftools sort {vcf} --temp-dir $BATCH_TMPDIR/ | bgzip -c > {bcftools_job.vcf_output['vcf.gz']};

    echo "indexing {bcftools_job.vcf_output['vcf.gz']}";
    tabix -p vcf {bcftools_job.vcf_output['vcf.gz']};
"""
    )
    b.write_output(
        bcftools_job.vcf_output, output_path(f'vcf/{version}/bctools/hail_filtered')
    )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
