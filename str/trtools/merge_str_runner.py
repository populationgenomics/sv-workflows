#!/usr/bin/env python3
"""
This script merges ExpansionHunter vcf.gz files into one combined VCF.
Please ensure merge_prep.py has been run on the vcf files prior to running mergeSTR.py

Optional ability to add in VCFs from another file directory (but must be sharded in the same way as the input-dir-1)

For example:
analysis-runner --access-level standard --dataset tob-wgs --description '5M merge TOB100' --output-dir 'str/5M_run_combined_vcfs/merge_str/v4' merge_str_runner.py --input-dir-1=gs://cpg-tob-wgs-main-analysis/str/5M_run_combined_vcfs/merge_str_prep/v4-2 --num-shards=50 --sample-list-1=gs://

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os

import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path


config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']


# inputs:


# input directory 1
@click.option('--input-dir-1', help='gs://...')
# input directory 2
@click.option('--input-dir-2', help='gs://...', default=None)
# sample list 1 (CPG sample IDs separated by \n)
@click.option('--sample-list-1', help='gs://...')
# sample list 2 (CPG sample IDs separated by \n)
@click.option('--sample-list-2', help='gs://...', default=None)
# input num shards
@click.option(
    '--num-shards',
    type=int,
    default=1,
    help=('Number of shards per sample to expect (if unsharded, choose 1)'),
)
@click.option(
    '--job-storage', help='Storage of the Hail batch job eg 30G', default='20G'
)
@click.command()
def main(
    job_storage, input_dir_1, input_dir_2, sample_list_1, sample_list_2, num_shards
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()
    for shard_index in enumerate(num_shards, 1):
        # Initialise TRTools job to run mergeSTR
        trtools_job = b.new_job(name=f'mergeSTR shard {shard_index}')
        trtools_job.image(TRTOOLS_IMAGE)
        trtools_job.cpu(8)
        trtools_job.storage(job_storage)
        trtools_job.declare_resource_group(
            vcf_output={
                'vcf': '{root}.vcf',
                'vcf.gz': '{root}.vcf.gz',
                'vcf.gz.tbi': '{root}.vcf.gz.tbi',
            }
        )

        # read in input file paths
        batch_vcfs = []
        num_samples = 0
        with to_path(sample_list_1).open() as f_1:
            ids_1 = [line.strip() for line in f_1]
            for id in ids_1:
                each_vcf = os.path.join(
                    input_dir_1, f'{id}_eh_shard{shard_index}.reheader.vcf.gz'
                )
                batch_vcfs.append(
                    b.read_input_group(
                        **{
                            'vcf.gz': each_vcf,
                            'vcf.gz.tbi': f'{each_vcf}.tbi',
                        }
                    )['vcf.gz']
                )
            num_samples = num_samples + len(ids_1)

        # if second file directory is specified, read in input file paths:
        if input_dir_2 is not None:
            with to_path(sample_list_2).open() as f_2:
                ids_2 = [line.strip() for line in f_2]
                for id in ids_2:
                    each_vcf = os.path.join(
                        input_dir_2, f'{id}_eh_shard{shard_index}.reheader.vcf.gz'
                    )
                    batch_vcfs.append(
                        b.read_input_group(
                            **{
                                'vcf.gz': each_vcf,
                                'vcf.gz.tbi': f'{each_vcf}.tbi',
                            }
                        )['vcf.gz']
                    )
                num_samples = num_samples + len(ids_2)

        trtools_job.command(
            f"""
        mergeSTR --vcfs {",".join(batch_vcfs)} --out {trtools_job.vcf_output} --vcftype eh
        bgzip -c {trtools_job.vcf_output}.vcf > {trtools_job.vcf_output['vcf.gz']}
        tabix -f -p vcf {trtools_job.vcf_output['vcf.gz']}  > {trtools_job.vcf_output['vcf.gz.tbi']}
        """
        )

        output_path_name = output_path(
            f'mergeSTR_{num_samples}_samples_eh_shard{shard_index}', 'analysis'
        )
        b.write_output(trtools_job.vcf_output['vcf.gz'], f'{output_path_name}.vcf.gz')
        b.write_output(
            trtools_job.vcf_output['vcf.gz.tbi'], f'{output_path_name}.vcf.gz.tbi'
        )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
