#!/usr/bin/env python3

"""
This script runs qcSTR() from TRTools package on a single/merged STR vcf file. 

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
# file-path
@click.option('--file-path', help='gs://...')
# caller
@click.option('--caller', help='gangstr or eh')

@click.command()
def main(
    file_path, caller
):  # pylint: disable=missing-function-docstring

    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    vcf_input = b.read_input(file_path)
    trtools_job = b.new_job(name="qcSTR")
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('20G')
    trtools_job.cpu(8)

    trtools_job.declare_resource_group(ofile = {'sample-callnum.pdf': '{root}-sample-callnum.pdf',
                                               'chrom-callnum.pdf': '{root}-chrom-callnum.pdf',
                                               'diffref-histogram.pdf': '{root}-diffref-histogram.pdf',
                                               'quality.pdf': '{root}-quality.pdf',
                                               'diffref-bias.pdf': '{root}-diffref-bias.pdf',
                                               'quality-per-locus.pdf': '{root}-quality-per-locus.pdf',
                                               'quality-sample-stratified.pdf': '{root}-quality-sample-stratified.pdf',
                                               'quality-per-sample.pdf': '{root}-quality-per-sample.pdf'
    })

    trtools_job.command(f"""
    set -ex;
    qcSTR --vcf {vcf_input} --vcftype {caller} --quality per-locus --quality sample-stratified --quality per-sample --out {trtools_job.ofile}
    
    """
    )
    num_samples = len(vcf_input)

    output_path_vcf = output_path(f'qCSTR_{num_samples}_samples_{caller}')
    b.write_output(trtools_job.ofile, output_path_vcf)

    b.run(wait = False)


if __name__ == '__main__':
    main()
