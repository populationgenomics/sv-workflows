#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script is step 4 of 5 for running associaTR.
It aims to:
- apply call level filters to the VCF (assumes dumpSTR locus level filters have already been applied)
specifically:
1) CI/CN <= 1 (CI = confidence interval size, CN = allele size)
2) Allele size - Mode allele size of the locus <=20 (filters out large expansions)

 analysis-runner --dataset "tob-wgs" \
    --description "Hail call level filters applied" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
    associatr_runner_hail.py --file-path=gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/dumpSTR/dumpSTR_filtered.vcf.gz

"""

import hail as hl
import click

from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch
from cpg_utils import to_path

from cpg_utils.hail_batch import output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']

def call_level_filter(file_path, ofile_path):
    hl.init()
    mt = hl.import_vcf(file_path,force_bgz = True)

    #wrangle CI (confidence interval columns)
    mt = mt.annotate_entries(allele_1_rep_length=hl.int(mt.REPCN.split("/")[0]))
    mt = mt.annotate_entries(allele_2_rep_length = hl.if_else(hl.len(mt.REPCN.split("/")) ==2, hl.int(mt.REPCN.split("/")[1]), hl.missing('int32')))
    mt = mt.annotate_entries(allele_1_REPCI=hl.str(mt.REPCI.split("/")[0]),
                            allele_2_REPCI=hl.if_else(hl.len(mt.REPCN.split("/")) ==2, hl.str(mt.REPCI.split("/")[1]), hl.missing('str')))
    mt = mt.annotate_entries(allele_1_REPCI_1=hl.int32(mt.allele_1_REPCI.split("-")[0]),
                            allele_1_REPCI_2=hl.int32(mt.allele_1_REPCI.split("-")[1]))
    mt = mt.annotate_entries(allele_2_REPCI_1 = hl.if_else(hl.len(mt.allele_2_REPCI.split("-")) ==2, hl.int(mt.allele_2_REPCI.split("-")[0]), hl.missing('int32')))
    mt = mt.annotate_entries(allele_2_REPCI_2 = hl.if_else(hl.len(mt.allele_2_REPCI.split("-")) ==2, hl.int(mt.allele_2_REPCI.split("-")[1]), hl.missing('int32')))
    mt = mt.annotate_entries(allele_1_REPCI_size=hl.int32(mt.allele_1_REPCI_2 - mt.allele_1_REPCI_1),
                            allele_2_REPCI_size=hl.int32(mt.allele_2_REPCI_2 - mt.allele_2_REPCI_1))
    mt = mt.annotate_entries(allele_1_REPCI_over_CN=hl.float64(mt.allele_1_REPCI_size/mt.allele_1_rep_length),
                            allele_2_REPCI_over_CN=hl.float64(mt.allele_2_REPCI_size/mt.allele_2_rep_length))
    mt = mt.drop('allele_1_REPCI_1', 'allele_1_REPCI_2', 'allele_2_REPCI_1', 'allele_2_REPCI_2','allele_1_REPCI', 'allele_2_REPCI','allele_1_REPCI_size', 'allele_2_REPCI_size')

    # Set all entry fields to missing where REP CI/CN >1
    mt = mt.annotate_entries(
        GT=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('call'), mt.GT),
        ADFL=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('str'), mt.ADFL),
        ADIR=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('str'), mt.ADIR),
        ADSP=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('str'), mt.ADSP),
        LC=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('float64'), mt.LC),
        REPCI=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('str'), mt.REPCI),
        REPCN=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('str'), mt.REPCN),
        SO=hl.if_else(mt.allele_1_REPCI_over_CN>1, hl.missing('str'), mt.SO),
        FILTER=hl.if_else(mt.allele_1_REPCI_over_CN>1, "NOCALL", mt.FILTER)
    )

    mt = mt.annotate_entries(
        GT=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('call'), mt.GT),
        ADFL=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('str'), mt.ADFL),
        ADIR=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('str'), mt.ADIR),
        ADSP=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('str'), mt.ADSP),
        LC=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('float64'), mt.LC),
        REPCI=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('str'), mt.REPCI),
        REPCN=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('str'), mt.REPCN),
        SO=hl.if_else(mt.allele_2_REPCI_over_CN>1, hl.missing('str'), mt.SO),
        FILTER=hl.if_else(mt.allele_2_REPCI_over_CN>1, "NOCALL", mt.FILTER)
    )

    # calculate mode allele per REPID (locus)
    mode_mt = mt.key_rows_by(REPID = mt.info.REPID)
    #annotate rows with counts of alleles that appear in each VARID
    ht = mode_mt.select_rows(
        alleles_rep_lengths = hl.agg.collect(mode_mt.allele_1_rep_length)
            .extend(hl.agg.collect(mode_mt.allele_2_rep_length))
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

    mt = mt.annotate_rows(mode = aggregated_table[mt.info.REPID])

    mt = mt.annotate_entries(allele_1_minus_mode = (mt.allele_1_rep_length- mt.mode.mode_allele),
                         allele_2_minus_mode = (mt.allele_2_rep_length-mt.mode.mode_allele))
    # Set all entry fields to missing where allele - mode allele > 20
    mt = mt.annotate_entries(
        GT=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('call'), mt.GT),
        ADFL=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('str'), mt.ADFL),
        ADIR=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('str'), mt.ADIR),
        ADSP=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('str'), mt.ADSP),
        LC=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('float64'), mt.LC),
        REPCI=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('str'), mt.REPCI),
        REPCN=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('str'), mt.REPCN),
        SO=hl.if_else(mt.allele_1_minus_mode>20, hl.missing('str'), mt.SO),
        FILTER=hl.if_else(mt.allele_1_minus_mode>20, "NOCALL", mt.FILTER)
    )

    mt = mt.annotate_entries(
        GT=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('call'), mt.GT),
        ADFL=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('str'), mt.ADFL),
        ADIR=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('str'), mt.ADIR),
        ADSP=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('str'), mt.ADSP),
        LC=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('float64'), mt.LC),
        REPCI=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('str'), mt.REPCI),
        REPCN=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('str'), mt.REPCN),
        SO=hl.if_else(mt.allele_2_minus_mode>20, hl.missing('str'), mt.SO),
        FILTER=hl.if_else(mt.allele_2_minus_mode>20, "NOCALL", mt.FILTER)
    )
    #drop unneccessary columns prior to writing out
    mt = mt.drop(mt.mode)
    mt = mt.drop(mt.low_depth_quality)
    mt = mt.drop('allele_1_rep_length','allele_2_rep_length','allele_1_REPCI_over_CN','allele_2_REPCI_over_CN','allele_1_minus_mode_imputed','allele_2_minus_mode_imputed')
    hl.export_vcf(mt, ofile_path)

@click.option(
    '--file-path',
    help='GCS file path to VCF (assumes dumpSTR locus level filters have already been applied)',
    type=str,
)
@click.command()
def main(file_path):

    b = get_batch()
    hail_job = b.new_python_job(name = f'Hail query call level filters')
    hail_job.image(config['workflow']['driver_image'])
    hail_job.storage('20G')
    hail_job.cpu(4)
    hail_job.call(call_level_filter, file_path, hail_job.ofile)

    bcftools_job = b.new_job(name = f'bgzip and tabix the Hail output VCF')
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.storage('20G')
    bcftools_job.cpu(4)

    bcftools_job.declare_resource_group(
        vcf_output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    bcftools_job.command(
        f"""
    set -ex;
    echo "Compressing";
    bcftools sort {hail_job.ofile} | bgzip -c > {bcftools_job.vcf_output['vcf.gz']};

    echo "indexing {bcftools_job.vcf_output['vcf.gz']}";
    tabix -p vcf {bcftools_job.vcf_output['vcf.gz']};
"""
    )
    b.write_output(bcftools_job.vcf_output, output_path(f'input_files/hail/hail_filtered'))

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter


