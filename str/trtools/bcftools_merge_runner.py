#!/usr/bin/env python3
# pylint: disable=duplicate-code,unsubscriptable-object
"""
This script merges ExpansionHunter VCFs together into one VCF using `bcftools merge`
Required input: --input-dir and internal sample CPG IDs
For example:
analysis-runner --access-level standard --dataset tob-wgs --description 'bcftools merge test' --output-dir 'str/5M_run_combined_vcfs/bcftools_merge/v4-3' bcftools_merge_runner.py --input-dir=gs://cpg-tob-wgs-main-analysis/str/5M_run_combined_vcfs/v4 CPGX CPGY

"""
import os
import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']

# inputs:
# input directory
@click.option('--input-dir', help='gs://...')
# input sample ID
@click.argument('internal-wgs-ids', nargs=-1)
@click.command()
def main(
    input_dir, internal_wgs_ids: list[str]
):
    """ Merge sample VCFs using bcftools merge """

    # Initializing Batch
    b = get_batch()

    # read in input file paths
    batch_vcfs = []
    for id in list(internal_wgs_ids):
        each_vcf = os.path.join(input_dir, f'{id}_eh.reheader.vcf.gz')
        batch_vcfs.append(b.read_input_group(**{
            'vcf.gz': each_vcf,
            'vcf.gz.tbi': f'{each_vcf}.tbi',
        })['vcf.gz'])
    num_samples = len(internal_wgs_ids)

    bcftools_job = b.new_job(name=f'Bcftools merge job')
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.cpu(4)
    bcftools_job.storage('20G')

    bcftools_job.declare_resource_group(vcf_out={'vcf.gz': '{root}.vcf.gz'})
    bcftools_job.command(
        f"""

        bcftools merge -m both -o {bcftools_job.vcf_out['vcf.gz']} -O z --threads 4 {" ".join(batch_vcfs)}

        """
    )
    # Output writing
    output_path_eh = output_path(f'merged_{num_samples}_eh', 'analysis')
    b.write_output(bcftools_job.vcf_sorted, output_path_eh)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,unsubscriptable-object
