#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses ExpansionHunterv5 to call STRs on WGS cram files.
Required input: --variant-catalog (file path to variant catalog), --dataset, and external sample IDs
For example:
analysis-runner --access-level test --dataset hgdp --description 'tester EH' --output-dir 'str/sensitivity-analysis/eh' eh_runner_hgdp.py --variant-catalog=gs://cpg-tob-wgs-test/hoptan-str/catalogs/sensitivity_analysis/marshfield_regions_eh.json --input-dir=gs://cpg-hgdp-test/cram/nagim

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
EH_IMAGE = config['images']['expansionhunter']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
@click.option('--input-dir', help='Full path to directory containing crams')
@click.command()
def main(
    variant_catalog, input_dir):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    #ref_fasta = ('gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta')
  
    ref_fasta = str(reference_path('broad/ref_fasta'))
    
    crams_path = []
    bucket_name, *components = input_dir[5:].split('/')

    client = storage.Client()

    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files: Set[str] = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    for file in files: 
        if file.endswith(".cram"): 
                crams_path.append(file)

    eh_regions = b.read_input(variant_catalog)

    # Iterate over each sample to call Expansion Hunter
    for cram_obj in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj, 'cram.crai': cram_obj+ '.crai'}
        )
        cpg_sample_id = cram_obj[30:38]

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

        # ExpansionHunter job initialisation
        eh_job = b.new_job(name=f'ExpansionHunter:{cpg_sample_id} running')
        eh_job.image(EH_IMAGE)
        eh_job.storage('50G')
        eh_job.cpu(8)

        eh_job.declare_resource_group(
            eh_output={
                'vcf': '{root}.vcf',
                'json': '{root}.json',
                'realigned_bam': '{root}_realigned.bam',
            }
        )

        eh_job.command(
            f"""
        ExpansionHunter  \\
        --reads {crams['cram']} \\
        --reference {ref.base} --variant-catalog {eh_regions}\\
         --threads 16 --analysis-mode streaming \\
         --output-prefix {eh_job.eh_output}
        """
        )
        # ExpansionHunter output writing
        eh_output_path = output_path(f'{cpg_sample_id}_eh')
        b.write_output(eh_job.eh_output, eh_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
