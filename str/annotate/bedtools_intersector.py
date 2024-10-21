#!/usr/bin/env python3



"""
Intersect with Arthur's annotations
analysis-runner --dataset "bioheart" --description "bedtools_intersector" --access-level "test" \
    --output-dir "str/arthur" bedtools_intersector.py \
    bedtools_intersector.py

"""
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

CATALOG_PATH = 'gs://cpg-bioheart-test/str/ncAnnot.v0.14.jul2024.bed'
TR_PATH = 'gs://cpg-bioheart-test/str/arthur/eTRs.bed'


# Initializing Batch
b = get_batch('Bedtools intersect')

# provision a new job
bedtools_job = b.new_job(name='Bedtools intersect')
bedtools_job.image(get_config()['images']['bedtools'])
bedtools_job.storage('20G')
bedtools_job.cpu(8)

# read input files
catalog = b.read_input(CATALOG_PATH)
tr = b.read_input(TR_PATH)

# set the job command
bedtools_job.command(f'bedtools intersect -a {tr} -b {catalog} -wo > {bedtools_job.ofile}')

# write output to GCP bucket for this dataset
b.write_output(bedtools_job.ofile, 'gs://cpg-bioheart-test/str/arthur/eTRs_intersect.bed')

b.run(wait=False)