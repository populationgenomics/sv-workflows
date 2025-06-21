#!/usr/bin/env python3



"""
Intersect with Arthur's annotations

analysis-runner --dataset "tenk10k" --description "bedtools_intersector" --access-level "test" \
    --output-dir "str/associatr/final_freeze/meta_fixed/arthur" bedtools_intersector.py \
    bedtools_intersector.py --tr-path=gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/arthur/input_files/fm_chr_pos.tsv

"""
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path
import click

CATALOG_PATH = 'gs://cpg-bioheart-test/str/ncAnnot.v0.14.jul2024.bed'

@click.option('--tr-path', type=str, help='Path to the tandem repeats file')
@click.command()
def main(tr_path):
    # Initializing Batch
    b = get_batch('Bedtools intersect')

    # provision a new job
    bedtools_job = b.new_job(name='Bedtools intersect')
    bedtools_job.image(get_config()['images']['bedtools'])
    bedtools_job.storage('20G')
    bedtools_job.memory('highmem')
    bedtools_job.cpu(16)

    # read input files
    catalog = b.read_input(CATALOG_PATH)
    tr = b.read_input(tr_path)
    tr_file_name = str(tr).split('/')[-1].split('.')[0]

    # set the job command
    bedtools_job.command(f'bedtools intersect -a {tr} -b {catalog} -wo > {bedtools_job.ofile}')

    # write output to GCP bucket for this dataset
    b.write_output(bedtools_job.ofile, output_path(f'o_files/{tr_file_name}_intersect.bed', 'analysis'))

    b.run(wait=False)

if __name__ == '__main__':
    main()