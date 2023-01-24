#!/usr/bin/env python3
# pylint: disable=import-error,duplicate-code

"""
This script uses GangSTRv2.5 to call STRs on WGS cram files.
Required input: --variant-catalog (file path to variant catalog), --dataset, and external sample IDs
For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'hoptan-str/CPG199760' gangstr_sample.py --variant-catalog=gs://cpg-tob-wgs-test/hoptan-str/catalogs/gangstr_pure_repeats_catalog.bed --dataset=tob-wgs

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

SAMTOOLS_IMAGE = config['images']['samtools']
GANGSTR_IMAGE = config['images']['gangstr']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
# input dataset
@click.option('--dataset', help='dataset to operate on, eg: tob-wgs')

@click.command()
def main(
    variant_catalog, dataset
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=config['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    if dataset == 'tob-wgs':
        ref_fasta = (
            'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta')
        
    crams_path = ["gs://cpg-tob-wgs-test/cram/nagim/CPG199760.cram"]

    gangstr_regions = b.read_input(variant_catalog)

    # Iterate over each sample to call GangSTR
    for cram_obj in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj['output'], 'cram.crai': cram_obj['output'] + '.crai'}
        )
        cpg_sample_id = "CPG199760"

        # Working with CRAM files requires the reference fasta
        ref = b.read_input_group(
            **dict(
                base=ref_fasta,
                fai=ref_fasta + '.fai',
                dict=ref_fasta.replace('.fasta', '')
                .replace('.fna', '')
                .replace('.fa', '')
                + '.dict',
            )
        )

        # GangSTR job initialisation
        gangstr_job = b.new_job(name=f'GangSTR:{cpg_sample_id} running')
        gangstr_job.image(GANGSTR_IMAGE)
        gangstr_job.storage('50G')
        gangstr_job.cpu(8)

        gangstr_job.declare_resource_group(
            gangstr_output={
                'vcf': '{root}.vcf',
                'insdata': '{root}.insdata.tab',
                'samplestats': '{root}.samplestats.tab',
            }
        )

        gangstr_job.command(
            f"""
        GangSTR --bam {crams['cram']} --ref {ref.base} --regions {gangstr_regions} --out {gangstr_job.gangstr_output} --bam-samps potato
        """
        )
        # GangSTR output writing
        gangstr_output_path = output_path(f'CPG199760_gangstr')
        b.write_output(gangstr_job.gangstr_output, gangstr_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
