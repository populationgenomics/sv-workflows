#!/usr/bin/env python3
"""
This script prepares GangSTR/EH VCF files for input into mergeSTR. 
Required input: --caller, --input-dir, and external sample IDs
For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester --output-dir 'tester' merge_prep.py --caller=eh --input-dir=gs://cpg-tob-wgs-main/str/expansionhunter/pure_repeats --dataset=tob-wgs TOBXXXX TOBXXXX

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

ref_fasta = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
BCFTOOLS_IMAGE = config['images']['bcftools']

# inputs:
# caller
@click.option('--caller', help='gangstr or eh')
# dataset
@click.option('--dataset', help='dataset eg tob-wgs')
# input directory
@click.option('--input-dir', help='gs://...')
# input sample ID
@click.argument('external-wgs-ids', nargs=-1)
@click.command()
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

    # Working with CRAM files requires the reference fasta
    ref = b.read_input_group(
        **dict(
            base=ref_fasta,
            fai=ref_fasta + '.fai',
            dict=ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )

    input_vcf_dict = {}

    if caller == 'eh':
        for id in list(external_id_to_cpg_id.values()):
            input_vcf_dict[id] = input_dir + "/" + id + "_EH.vcf"
    elif caller == 'gangstr':
        for id in list(external_id_to_cpg_id.values()):
            input_vcf_dict[id] = input_dir + "/" + id + "_gangstr.vcf"
    else:
        raise Exception("Invalid caller")

    for id in list(input_vcf_dict.keys()):

        bcftools_job = b.new_job(name=f'{id} {caller} Files prep')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.storage('20G')
        bcftools_job.cpu(8)

        vcf_input = b.read_input(input_vcf_dict[id])

        if caller == "eh":
            bcftools_job.declare_resource_group(
                vcf_sorted={
                    'vcf.gz': '{root}.vcf.gz',
                    'reheader.vcf.gz': '{root}.reheader.vcf.gz',
                    'vcf.gz.tbi': '{root}.reheader.vcf.gz.tbi',
                }
            )
            bcftools_job.command(
                f"""

                bgzip -c {vcf_input} > {bcftools_job.vcf_sorted['vcf.gz']}
            
                bcftools reheader -f {ref.fai} -o {bcftools_job.vcf_sorted['reheader.vcf.gz']} {bcftools_job.vcf_sorted['vcf.gz']} 

                tabix -f -p vcf {bcftools_job.vcf_sorted['reheader.vcf.gz']} > {bcftools_job.vcf_sorted['vcf.gz.tbi']}
            
                """
            )
            # Output writing
            output_path_eh = output_path(f'{id}_eh')
            b.write_output(bcftools_job.vcf_sorted['reheader.vcf.gz'], output_path_eh+".reheader.vcf.gz")
            b.write_output(bcftools_job.vcf_sorted['vcf.gz.tbi'], output_path_eh+".vcf.gz.tbi")

        else:
            bcftools_job.declare_resource_group(
                vcf_sorted={
                    'vcf.gz': '{root}.vcf.gz',
                    'vcf.gz.tbi': '{root}.reheader.vcf.gz.tbi',
                }
            )
            bcftools_job.command(
                f"""

                bcftools sort {vcf_input} | bgzip -c  > {bcftools_job.vcf_sorted['vcf.gz']}
            
                tabix -f -p vcf {bcftools_job.vcf_sorted['vcf.gz']} > {bcftools_job.vcf_sorted['vcf.gz.tbi']}
            
                """
            )
            # Output writing
            output_path_gangstr = output_path(f'{id}_gangstr')
            b.write_output(bcftools_job.vcf_sorted['vcf.gz'], output_path_gangstr+".vcf.gz")
            b.write_output(bcftools_job.vcf_sorted['vcf.gz.tbi'], output_path_gangstr+".vcf.gz.tbi")

    b.run(wait=False)


if __name__ == '__main__':
    main()
