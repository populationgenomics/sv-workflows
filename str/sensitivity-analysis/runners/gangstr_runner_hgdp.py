#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses GangSTR to call STRs on WGS cram files.
For example:
analysis-runner --access-level test --dataset hgdp --description 'gangstr run' --output-dir 'str/sensitivity-analysis/gangstr/untrimmed_coordinates' gangstr_runner_hgdp.py --variant-catalog=gs://cpg-hgdp-test/str/untrimmed_coordinates_resources/marshfield_regions_gangstr_untrimmed_coordinates.bed 

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import logging

import click
from sample_metadata.apis import SampleApi, ParticipantApi

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, reference_path
from cpg_workflows.batch import get_batch
from google.cloud import storage

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
GANGSTR_IMAGE = config['images']['gangstr']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
@click.command()
def main(
    variant_catalog):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    ref_fasta = ('gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta')
  
    #ref_fasta = str(reference_path('broad/ref_fasta'))
    
    crams_path = []
    """
    bucket_name, *components = input_dir[5:].split('/')

    client = storage.Client()

    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files: Set[str] = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    for file in files: 
        if file.endswith(".cram"): 
                crams_path.append(file)
    """

    data = ParticipantApi().get_participants('hgdp-test')
    population_groups = ["French", "Yoruba", "Han", "Northern Han"]
    participant_ids = []

    for i in data: 
        if i["meta"]["Population name"] in population_groups: 
            participant_ids.append(i['id']) 
    samples = SampleApi().get_samples(body_get_samples={'project_ids':['hgdp-test'], "participant_ids": participant_ids})
    for i in samples: 
        crams_path.append("gs://cpg-hgdp-test/cram/"+i["id"]+".cram")

    gangstr_regions = b.read_input(variant_catalog)

    # Iterate over each sample to call Expansion Hunter
    for cram_obj in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj, 'cram.crai': cram_obj+ '.crai'}
        )
        cpg_sample_id = cram_obj.replace('.cram','')[24:]

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
        GangSTR --bam {crams['cram']} --ref {ref.base} --regions {gangstr_regions} --out {gangstr_job.gangstr_output} --bam-samps {cpg_sample_id}
        """
        )
        # GangSTR output writing
        gangstr_output_path = output_path(f'{cpg_sample_id}_gangstr')
        b.write_output(gangstr_job.gangstr_output, gangstr_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter