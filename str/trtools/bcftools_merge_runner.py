#!/usr/bin/env python3
# pylint: disable=duplicate-code
"""
This script merges ExpansionHunter VCFs together into one VCF using `bcftools merge`
Required input: --input-dir and external sample IDs
For example:
analysis-runner --access-level standard --dataset tob-wgs --description 'bcftools merge test ' --output-dir 'str/5M_run_combined_vcfs/bcftools_merge/v4-3' bcftools_merge_runner.py --input-dir=gs://cpg-tob-wgs-main-analysis/str/5M_run_combined_vcfs/v4 CPGtestersorted2

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path, reference_path

config = get_config()

REF_FASTA = str(reference_path('broad/ref_fasta'))
BCFTOOLS_IMAGE = config['images']['bcftools']


# inputs:
# input directory
@click.option('--input-dir', help='gs://...')
# input sample ID
@click.argument('internal-wgs-ids', nargs=-1)
@click.command()
def main(
    input_dir, internal_wgs_ids: list[str]
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    # Working with CRAM files requires the reference fasta
    ref = b.read_input_group(
        **dict(
            base=REF_FASTA,
            fai=REF_FASTA + '.fai',
            dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )

    # read in input file paths
    vcffuse_path = []
    for id in list(internal_wgs_ids):
        vcf = os.path.join(input_dir, f'{id}_eh.reheader.vcf.gz')
        suffix = vcf.removeprefix('gs://').split('/', maxsplit=1)[1]
        vcffuse_path.append(f'/vcffuse/{suffix}')
    num_samples = len(vcffuse_path)
    vcffuse_path = ' '.join(vcffuse_path)  # string format for input into bcftools merge


    bcftools_job = b.new_job(name=f'Bcftools merge job')
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.cpu(4)
    bcftools_job.storage('15G')


    bcftools_job.declare_resource_group(
        vcf_out={
            'vcf.gz': '{root}.vcf.gz'
                            }
    )
    bcftools_job.command(
        f"""

        bcftools merge -m both -o {bcftools_job.vcf_out} -O z --threads 4 {vcffuse_path}

        """
    )
    # Output writing
    output_path_eh = output_path(f'{id}_eh', 'analysis')
    b.write_output(bcftools_job.vcf_sorted, output_path_eh)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
