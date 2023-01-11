#!/usr/bin/env python3

"""
Converts BED file to FASTA sequence 
"""
import os
import hailtop.batch as hb
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import (
    remote_tmpdir,
)
config = get_config()
DATASET = os.getenv('DATASET')
HAIL_BUCKET = os.getenv('HAIL_BUCKET')
OUTPUT_SUFFIX = os.getenv('OUTPUT')
BILLING_PROJECT = os.getenv('HAIL_BILLING_PROJECT')
ACCESS_LEVEL = os.getenv('ACCESS_LEVEL')

CATALOG_PATH = 'gs://cpg-tob-wgs-test/hoptan-str/bed_catalog_without_complex_repeats.bed.txt'
BEDTOOLS_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bedtools:v2.30.0'
REF_FASTA = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'

@click.command()

def main():  # pylint: disable=missing-function-docstring 

    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))   
    bedtools_job = b.new_job(name = f'Get sequence data for FASTA files')
    catalog = b.read_input(CATALOG_PATH)
    fasta = b.read_input(REF_FASTA)
    bedtools_job.image(BEDTOOLS_IMAGE)
    bedtools_job.storage('20G')
    bedtools_job.cpu(8)

    bedtools_job.command(f"""

        bedtools getfasta -fi {fasta} -bed {catalog} > {bedtools_job.ofile}
        
        """)

    #Output writing 
    out_fname = f'catalog_fasta_sequences.fasta'
    output_path = f'gs://cpg-tob-wgs-test/hoptan-str/{out_fname}'
    b.write_output(bedtools_job.ofile, output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
