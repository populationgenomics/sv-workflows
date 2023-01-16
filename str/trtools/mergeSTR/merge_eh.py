#!/usr/bin/env python

"""
This script merges GangSTR or ExpansionHunter vcf.gz files into one combined VCF. 
Please ensure merge_prep.py has been run on the vcf files prior to running mergeSTR.py

For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester --output-dir 'tester' mergeSTR.py --caller=eh --input-dir=gs://cpg-tob-wgs-main/str/expansionhunter/pure_repeats --dataset=tob-wgs TOBXXXX TOBXXXX

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os
import logging

import click
import hailtop.batch as hb

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.apis import AnalysisApi, SampleApi
from sample_metadata.models import AnalysisStatus

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path, reference_path

config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']

# inputs:
# caller
@click.option('--caller', help='gangstr or eh')
# dataset
@click.option('--dataset', help='dataset eg tob-wgs')
# input directory
@click.option('--input-dir', help='gs://...')
# input sample ID
@click.argument('external-wgs-ids', nargs=-1)

def main(
    dataset, caller, input_dir, external_wgs_ids: list[str]
):  # pylint: disable=missing-function-docstring

    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    
    external_id_to_cpg_id: dict[str, str] = SampleApi().get_sample_id_map_by_external(
        dataset, list(external_wgs_ids)
    )
    cpg_id_to_external_id = {
        cpg_id: external_wgs_id
        for external_wgs_id, cpg_id in external_id_to_cpg_id.items()
    }

    input_file = b.read_input_group(
    TOB01784 = "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB01784_EHreheader.vcf.gz",
    TOB01784_tbi = "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB01784_EHreheader.vcf.gz.tbi",
    TOB01791= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB01791_EHreheader.vcf.gz",
    TOB01791_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB01791_EHreheader.vcf.gz.tbi",
    TOB0901= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB0901_EHreheader.vcf.gz",
    TOB0901_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB0901_EHreheader.vcf.gz.tbi",
    TOB1086= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1086_EHreheader.vcf.gz",
    TOB1086_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1086_EHreheader.vcf.gz.tbi",
    TOB1170= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1170_EHreheader.vcf.gz",
    TOB1170_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1170_EHreheader.vcf.gz.tbi",
    TOB1213= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1213_EHreheader.vcf.gz",
    TOB1213_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1213_EHreheader.vcf.gz.tbi",
    TOB1458= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1458_EHreheader.vcf.gz",
    TOB1458_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1458_EHreheader.vcf.gz.tbi",
    TOB1567= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1567_EHreheader.vcf.gz",
    TOB1567_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1567_EHreheader.vcf.gz.tbi",
    TOB1578= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1578_EHreheader.vcf.gz",
    TOB1578_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1578_EHreheader.vcf.gz.tbi",
    TOB1751= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1751_EHreheader.vcf.gz",
    TOB1751_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/TOB1751_EHreheader.vcf.gz.tbi"
)

    trtools_job = b.new_job(name = "mergeSTR")
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('200G')
    trtools_job.cpu(16)

    trtools_job.declare_resource_group(ofile = {'vcf': '{root}.vcf'})
        
    trtools_job.command(f"""
     
    mergeSTR --vcfs {input_file.TOB01784},{input_file.TOB01791},{input_file.TOB0901},{input_file.TOB1086},{input_file.TOB1170},{input_file.TOB1213},{input_file.TOB1458},{input_file.TOB1567},{input_file.TOB1578},{input_file.TOB1751} --out {trtools_job.ofile} --vcftype eh
     
    """)

    eh_out_name = "tob10_EH_merge"
    eh_output_path = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/merge/{eh_out_name}'
    b.write_output(trtools_job.ofile, eh_output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
