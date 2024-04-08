#!/usr/bin/env python3

"""
Converts coordinates in BED file to FASTA sequence
"""
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path


CATALOG_PATH = (
    'gs://cpg-tob-wgs-test/hoptan-str/bed_catalog_without_complex_repeats.bed.txt'
)
REF_FASTA = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'


# Initializing Batch
b = get_batch('FASTA conversion')

# provision a new job
bedtools_job = b.new_job(name=f'Get sequence data for FASTA files')
bedtools_job.image(get_config()['images']['bedtools'])
bedtools_job.storage('20G')
bedtools_job.cpu(8)

# read input files
catalog = b.read_input(CATALOG_PATH)
fasta = b.read_input(REF_FASTA)

# set the job command
bedtools_job.command(
    f'bedtools getfasta -fi {fasta} -bed {catalog} > {bedtools_job.ofile}'
)

# write output to GCP bucket for this dataset
b.write_output(
    bedtools_job.ofile, output_path('Illumina_catalog_sequences.fasta.txt', 'analysis')
)

b.run(wait=False)
