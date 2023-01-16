#!/usr/bin/env python3

"""
This script runs compareSTR() from TRTools package to compare calls between 2 STR VCFs.

For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' compareSTR_runner.py --file-path-1=gs://cpg-tob-wgs-test/hoptan-str/mergeSTR/mergeSTR_2_samples_gangstr.vcf
--file-path-2=gs://cpg-tob-wgs-test/hoptan-str/mergeSTR/mergeSTR_2_samples_eh.vcf --caller-1=gangstr --caller-2=eh

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os
import logging

import click
import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']

# inputs:
# file-path-1
@click.option('--file-path-1', help='gs://...to VCF 1')
# file-path-2
@click.option('--file-path-2', help='gs://... to VCF 2')
# caller-1
@click.option('--caller-1', help='gangstr or eh')
# caller-2
@click.option('--caller-2', help='gangstr or eh')
@click.command()

def main(file_path_1, file_path_2, caller_1, caller_2):  # pylint: disable=missing-function-docstring

    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    vcf_input_1 = b.read_input(file_path_1)
    vcf_input_2 = b.read_input(file_path_2)
    trtools_job = b.new_job(name=f"compareSTR")
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('20G')
    trtools_job.cpu(8)

    trtools_job.declare_resource_group(ofile = {'overall.tab': '{root}-overall.tab',
                                               'bubble-periodALL.pdf': '{root}-bubble-periodALL.pdf',
                                               'locuscompare.tab': '{root}-locuscompare.tab',
                                               'locuscompare.pdf': '{root}-locuscompare.pdf',
                                               'samplecompare.tab': '{root}-samplecompare.tab',
                                               'samplecompare.pdf': '{root}-samplecompare.pdf'


    })

    trtools_job.command(f"""
    set -ex;
    compareSTR --vcf1 {vcf_input_1} --vcf2 {vcf_input_2} --vcftype1 {caller_1} --vcftype2 {caller_2} --out {trtools_job.ofile}
    
    """
    )
    output_path_vcf = output_path(f'compareSTR_samples_{caller_1}_{caller_2}')
    b.write_output(trtools_job.ofile, output_path_vcf)

    b.run(wait=False)


if __name__ == '__main__':
    main()

