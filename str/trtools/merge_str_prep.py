#!/usr/bin/env python3
# pylint: disable=duplicate-code
"""
This script prepares GangSTR/EH VCF files for input into mergeSTR. 
Required input: --caller, --input-dir, and external sample IDs
For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester --output-dir 'tester' merge_prep.py --caller=eh --input-dir=gs://cpg-tob-wgs-main/str/expansionhunter/pure_repeats --dataset=tob-wgs CPGXXXX CPGXXXX

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch


config = get_config()

REF_FASTA = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
BCFTOOLS_IMAGE = config['images']['bcftools']


# inputs:
# caller
@click.option(
    '--caller',
    help='gangstr or eh',
    type=click.Choice(['eh', 'gangstr'], case_sensitive=True),
)
# input directory
@click.option('--input-dir', help='gs://...')
# input sample ID
@click.argument('internal-wgs-ids', nargs=-1)
@click.command()
def main(
    caller, input_dir, internal_wgs_ids: list[str]
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

    input_vcf_dict = {}

    input_vcf_dict = {id: os.path.join(input_dir, f'{id}_{caller}.vcf') for id in internal_wgs_ids}

    for id in list(input_vcf_dict.keys()):
        bcftools_job = b.new_job(name=f'{id} {caller} Files prep')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.storage('20G')
        bcftools_job.cpu(8)

        vcf_input = b.read_input(input_vcf_dict[id])

        if caller == 'eh':
            bcftools_job.declare_resource_group(
                vcf_sorted={
                    'vcf.gz': '{root}.vcf.gz',
                    'reheader.vcf.gz': '{root}.reheader.vcf.gz',
                    'reheader.vcf.gz.tbi': '{root}.reheader.vcf.gz.tbi',
                }
            )
            bcftools_job.command(
                f"""

                bgzip -c {vcf_input} > {bcftools_job.vcf_sorted['vcf.gz']}
            
                bcftools reheader -f {ref.fai} -o {bcftools_job.vcf_sorted['reheader.vcf.gz']} {bcftools_job.vcf_sorted['vcf.gz']} 

                tabix -f -p vcf {bcftools_job.vcf_sorted['reheader.vcf.gz']} 
            
                """
            )
            # Output writing
            output_path_eh = output_path(f'{id}_eh', 'analysis')
            b.write_output(bcftools_job.vcf_sorted, output_path_eh)

        else:
            bcftools_job.declare_resource_group(
                vcf_sorted={
                    'vcf.gz': '{root}.vcf.gz',
                    'vcf.gz.tbi': '{root}.vcf.gz.tbi',
                }
            )
            bcftools_job.command(
                f"""

                bcftools sort {vcf_input} | bgzip -c  > {bcftools_job.vcf_sorted['vcf.gz']}
            
                tabix -p vcf {bcftools_job.vcf_sorted['vcf.gz']}
            
                """
            )
            # Output writing
            output_path_gangstr = output_path(f'{id}_gangstr', 'analysis')
            b.write_output(bcftools_job.vcf_sorted, output_path_gangstr)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
