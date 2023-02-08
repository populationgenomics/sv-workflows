#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses GangSTR to call STRs on WGS cram files.
For example:
analysis-runner --access-level test --dataset hgdp --description 'hipstr run' --output-dir 'str/sensitivity-analysis/hipstr/untrimmed_coordinates' hipstr_runner_hgdp.py --variant-catalog=gs://cpg-hgdp-test/str/untrimmed_coordinates_resources/hg38_hipstr_catalog_untrimmed_coordinates.bed --input-dir=gs://cpg-hgdp-test/cram/nagim

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import logging

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, reference_path
from cpg_workflows.batch import get_batch
from google.cloud import storage

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
HIPSTR_IMAGE = config['images']['hipstr']
BCFTOOLS_IMAGE = config['images']['bcftools']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
@click.option('--input-dir', help='Full path to directory containing crams')
@click.command()
def main(
    variant_catalog, input_dir):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    ref_fasta = ('gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta')
  
    #ref_fasta = str(reference_path('broad/ref_fasta'))
    
    crams_path = []
    bucket_name, *components = input_dir[5:].split('/')

    client = storage.Client()

    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files: Set[str] = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    for file in files: 
        if file.endswith(".cram"): 
                crams_path.append(file)

    gangstr_regions = b.read_input(variant_catalog)

    # Iterate over each sample to call Expansion Hunter
    for cram_obj in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj, 'cram.crai': cram_obj+ '.crai'}
        )
        cpg_sample_id = cram_obj.replace('.cram','')[30:]

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

        hipstr_job = b.new_job(name=f'HipSTR:{cpg_sample_id} running')
        hipstr_job.image(HIPSTR_IMAGE)
        hipstr_job.storage('50G')
        hipstr_job.cpu(8)

        hipstr_job.declare_resource_group(
            hipstr_output={
                'vcf.gz': '{root}.vcf.gz',
                'viz.gz': '{root}.viz.gz'
            }
        )

        hipstr_job.command(
            f"""
        HipSTR --bams {crams['cram']} --fasta {ref.base} --regions {gangstr_regions} --str-vcf {hipstr_job.hipstr_output['vcf.gz']} --viz-out {hipstr_job.hipstr_output['viz.gz']} --min-reads 25
        """
        )
        # HipSTR output writing
        hipstr_output_path_vcf = output_path(f'{cpg_sample_id}_hipstr.vcf.gz')
        b.write_output(hipstr_job.hipstr_output['vcf.gz'], hipstr_output_path_vcf)

        hipstr_output_path_viz = output_path(f'{cpg_sample_id}_hipstr.viz.gz')
        b.write_output(hipstr_job.hipstr_output['viz.gz'], hipstr_output_path_viz)

        bcftools_job = b.new_job(name=f'{cpg_sample_id} Unzip VCF')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.storage('20G')
        bcftools_job.cpu(8)
        bcftools_job.command(
                f"""
                bgzip -d {hipstr_job.hipstr_output['vcf.gz']} > {bcftools_job.ofile}
            
                """
            )
        bcftools_job_output_path = output_path(f'{cpg_sample_id}_hipstr.vcf')
        b.write_output(bcftools_job.ofile, bcftools_job_output_path)


    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
