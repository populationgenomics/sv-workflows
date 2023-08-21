#!/usr/bin/env python3
"""
This script merges ExpansionHunter vcf.gz files into one combined VCF. 
Please ensure merge_prep.py has been run on the vcf files prior to running mergeSTR.py

For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'str/expansionhunter/test_chr_22/mergeSTR' merge_str_runner.py --input-dir=gs://cpg-tob-wgs-test-analysis/str/expansionhunter/test_chr_22/mergeSTRprep --dataset=tob-wgs CPG199760 CPG199778

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os

import click

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
@click.argument('internal-wgs-ids', nargs=-1)
@click.command()
def main(
    dataset, input_dir, internal_wgs_ids: list[str]
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    # Initialise TRTools job to run mergeSTR
    trtools_job = b.new_job(name='mergeSTR')
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.cpu(8)
    # mount using cloudfuse for reading input files
    trtools_job.cloudfuse(f'cpg-{dataset}-main-analysis', '/vcffuse')
    trtools_job.declare_resource_group(vcf_output={'vcf': '{root}.vcf','vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})


    # read in input file paths
    vcffuse_path = []
    for id in list(internal_wgs_ids):
        vcf = os.path.join(input_dir, f'{id}_eh.reheader.vcf.gz')
        suffix = vcf.removeprefix('gs://').split('/', maxsplit=1)[1]
        vcffuse_path.append(f'/vcffuse/{suffix}')
    num_samples = len(vcffuse_path)
    vcffuse_path = ','.join(vcffuse_path)  # string format for input into mergeSTR

    trtools_job.command(
        f"""
    mergeSTR --vcfs {vcffuse_path} --out {trtools_job.vcf_output} --vcftype eh
    bgzip -c {trtools_job.vcf_output}.vcf > {trtools_job.vcf_output['vcf.gz']}
    tabix -p vcf {trtools_job.vcf_output['vcf.gz']}
    """
    )

    output_path_name = output_path(f'mergeSTR_{num_samples}_samples_eh', 'analysis')
    b.write_output(trtools_job.vcf_output['vcf.gz'], output_path_name)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
