#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals
"""
This script merges two (or more) mergeSTR outputs into a new mergeSTR output (rolling merge).
Specify as args the paths to the mergeSTR output (including file name prefix).

For example:
analysis-runner --access-level test --dataset bioheart --description '5M-3M mergeSTR tester' --output-dir 'str/polymorphic_run/merge_str/v1_n2412' rolling_merge_str_runner.py --num-shards=27 \
gs://cpg-bioheart-test-analysis/str/polymorphic_run/merge_str/bioheart/v2_n367/mergeSTR_367_samples_eh gs://cpg-bioheart-test-analysis/str/polymorphic_run_n2045/merge_str/v1/mergeSTR_2045_samples_eh

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
# num shards
@click.option(
    '--num-shards',
    type=int,
    default=1,
    help=('Number of shards per sample to expect (if unsharded, choose 1)'),
)
@click.option('--job-storage', help='Storage of the Hail batch job eg 30G', default='20G')
@click.option('--job-memory', help='Memory of the Hail batch job', default='standard')
@click.option('--job-cpu', help='Number of CPUs of the Hail batch job', default=8)
# input sample ID
@click.argument('input-file-paths', nargs=-1)
@click.command()
def main(
    job_storage,
    job_memory,
    job_cpu,
    num_shards,
    input_file_paths,
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()
    for shard_index in range(1, num_shards + 1):
        # Initialise TRTools job to run mergeSTR
        trtools_job = b.new_job(name=f'mergeSTR shard {shard_index}')
        trtools_job.image(TRTOOLS_IMAGE)
        trtools_job.cpu(job_cpu)
        trtools_job.memory(job_memory)
        trtools_job.storage(job_storage)
        trtools_job.declare_resource_group(
            vcf_output={
                'vcf': '{root}.vcf',
                'vcf.gz': '{root}.vcf.gz',
                'vcf.gz.tbi': '{root}.vcf.gz.tbi',
            },
        )

        # read in input file paths
        batch_vcfs = []
        num_samples = 0
        input_file_paths = list(input_file_paths)
        for input_file_path in input_file_paths:
            each_vcf = input_file_path + f'_shard{shard_index}.vcf.gz'
            batch_vcfs.append(
                b.read_input_group(
                    **{
                        'vcf.gz': each_vcf,
                        'vcf.gz.tbi': f'{each_vcf}.tbi',
                    },
                )['vcf.gz'],
            )
        num_samples = sum(int(path.split('_')[-3]) for path in input_file_paths)

        batch_vcfs_list = '${BATCH_TMPDIR}/batch_vcfs.list'

        newline = '\n'
        # Writing cat...EOF on one line avoids any indentation on the batch_vcfs and EOF lines
        trtools_job.command(
            f"""
        cat <<EOF >{batch_vcfs_list}\n{newline.join(batch_vcfs)}\nEOF

        mergeSTR --vcfs-list {batch_vcfs_list} --out {trtools_job.vcf_output} --vcftype eh
        bgzip -c {trtools_job.vcf_output}.vcf > {trtools_job.vcf_output['vcf.gz']}
        tabix -f -p vcf {trtools_job.vcf_output['vcf.gz']}  > {trtools_job.vcf_output['vcf.gz.tbi']}
        """,
        )

        output_path_name = output_path(f'mergeSTR_{num_samples}_samples_eh_shard{shard_index}', 'analysis')
        b.write_output(trtools_job.vcf_output['vcf.gz'], f'{output_path_name}.vcf.gz')
        b.write_output(trtools_job.vcf_output['vcf.gz.tbi'], f'{output_path_name}.vcf.gz.tbi')

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
