#!/usr/bin/env python3
# pylint: disable=duplicate-code
"""
This script runs statSTR() from TRTools package on a single/merged STR vcf file and outputs various statistics

For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' statSTR_runner.py --caller=eh --file-path=gs://cpg-tob-wgs-test/hoptan-str/mergeSTR/mergeSTR_2_samples_gangstr.vcf

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os

import click
import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']


# inputs:
# file-path
@click.option('--file-path', help='gs://...')
# caller
@click.option('--caller', help='gangstr or eh')
@click.command()
def main(file_path, caller):  # pylint: disable=missing-function-docstring

    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    vcf_input = b.read_input(file_path)
    trtools_job = b.new_job(name=f'statSTR {caller}')

    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('20G')
    trtools_job.cpu(8)

    trtools_job.declare_resource_group(ofile={'tab': '{root}.tab'})

    trtools_job.command(
        f"""
        set -ex;
        statSTR --vcf {vcf_input} --vcftype {caller} --out {trtools_job.ofile} --thresh --afreq --acount --hwep --het --entropy --mean --mode --var --numcalled
    
        """
    )

    output_path_vcf = output_path(f'statSTR_samples_{caller}')
    b.write_output(trtools_job.ofile, output_path_vcf)

    b.run(wait=False)


if __name__ == '__main__':
    main()# pylint: disable=no-value-for-parameter
