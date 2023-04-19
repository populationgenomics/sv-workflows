#!/usr/bin/env python3
"""
This script merges ExpansionHunter vcf.gz files into one combined VCF. 
Please ensure merge_prep.py has been run on the vcf files prior to running mergeSTR.py

For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester --output-dir 'tester' merge_str_runner.py --input-dir=gs://cpg-tob-wgs-main/str/expansionhunter/pure_repeats --dataset=tob-wgs TOBXXXX TOBXXXX

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os

import click

from sample_metadata.apis import SampleApi

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch


config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']


# inputs:

# dataset
@click.option('--dataset', help='dataset eg tob-wgs')
# input directory
@click.option('--input-dir', help='gs://...')
# input sample ID
@click.argument('external-wgs-ids', nargs=-1)
@click.command()
def main(
    dataset, input_dir, external_wgs_ids: list[str]
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    external_id_to_cpg_id: dict[str, str] = SampleApi().get_sample_id_map_by_external(
        dataset, list(external_wgs_ids)
    )

    trtools_job = b.new_job(name='mergeSTR')
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('20G')
    trtools_job.cpu(8)
    trtools_job.cloudfuse(f'cpg-{dataset}-main', '/vcffuse')

    vcffuse_path = []
    for id in list(external_id_to_cpg_id.values()):
        vcf = os.path.join(input_dir, f'{id}_eh.reheader.vcf.gz')
        suffix = vcf.removeprefix('gs://').split('/', maxsplit=1)[1]
        vcffuse_path.append(f'/vcffuse/{suffix}')
    vcffuse_path = ','.join(vcffuse_path)  # string format for input into mergeSTR

    trtools_job.declare_resource_group(ofile={'vcf': '{root}.vcf'})

    trtools_job.command(
        f"""
     
    mergeSTR --vcfs {vcffuse_path} --out {trtools_job.ofile} --vcftype eh
     
    """
    )
    num_samples = len(vcffuse_path)

    output_path_vcf = output_path(f'mergeSTR_{num_samples}_samples_eh')
    b.write_output(trtools_job.ofile, output_path_vcf)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
