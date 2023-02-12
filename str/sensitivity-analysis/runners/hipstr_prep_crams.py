#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script subsets crams for regions that intersect the HipSTR variant catalog +/- 10,000bp
For example:
analysis-runner --access-level test --dataset hgdp --description 'hipstr subset cram' --output-dir 'str/sensitivity-analysis/hipstr/subset_crams' hipstr_prep_crams.py --variant-catalog=gs://cpg-hgdp-test/str/untrimmed_coordinates_resources/hg38_hipstr_catalog_untrimmed_coordinates_1_based_10000bp_padding.bed

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import logging
import hailtop.batch as hb

import click
from sample_metadata.apis import SampleApi, ParticipantApi
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
click.command()
def main(
    variant_catalog):  # pylint: disable=missing-function-docstring
    # Initializing Batch 
    b = get_batch()

    ref_fasta = ('gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta')
  
    #ref_fasta = str(reference_path('broad/ref_fasta'))
    crams_gcp_path =["gs://cpg-hgdp-test/cram/nagim/CPG198697.cram"]
    
    hipstr_regions = b.read_input(variant_catalog)

    cram_collection = {}

    # Iterate over each sample to add to cram_collection
    for cram_obj in crams_gcp_path:
        cpg_sample_id = cram_obj.replace('.cram','')[30:]
        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj, 'cram.crai': cram_obj + '.crai'}
        )
        cram_collection[cpg_sample_id] = crams['cram']
    
    crams_batch_path = ""
    for i in cram_collection:
        crams_batch_path+=str(cram_collection[i])
        crams_batch_path+= ","
    crams_batch_path = crams_batch_path[:-1]

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
    samtools_job = b.new_job(name=f'Samtools running')
    samtools_job.image(SAMTOOLS_IMAGE)
    samtools_job.storage('50G')
    samtools_job.cpu(8)

    samtools_job.command(
        f"""
    samtools view -C -T {ref.base} --region-file {hipstr_regions} > {samtools_job.ofile}
    """
    )
    samtools_job_output_path = output_path(f'tester_hipstr_subset.cram')
    b.write_output(samtools_job.ofile, samtools_job_output_path)


    b.run(wait = False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
