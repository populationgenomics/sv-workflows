#!/usr/bin/env python3
# pylint: disable=duplicate-code
"""
This script runs statSTR() from TRTools package on a single/merged STR vcf file and outputs various statistics
FYI there is a bug in statSTR preventing it from aggregating stats for chrY loci so please filter out chrY loci prior.

For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'hoptan-str/shard_workflow_test/sharded_stat_str' stat_str_runner.py --caller=eh --input-dir=gs://cpg-tob-wgs-test/hoptan-str/shard_workflow_test/merge_str --sharded

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""

import click

from cpg_utils.config import get_config
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path


config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']


# inputs:
# file-path
@click.option(
    '--input-dir', help='gs://...file path if unsharded, input directory if sharded'
)
# caller
@click.option(
    '--caller',
    help='gangstr or eh',
    type=click.Choice(['eh', 'gangstr'], case_sensitive=True),
)
@click.option(
    '--job-storage', help='Storage of the Hail batch job eg 30G', default='16G'
)
# sharded flag
@click.option('--sharded', is_flag=True, help='Assume a sharded catalog was used')
@click.option('--job-memory', help='Memory of the Hail batch job eg 64G', default='16G')
@click.command()
def main(
    input_dir, caller, job_storage, job_memory, sharded
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    if sharded:
        input_vcf_list = [
            str(gspath) for gspath in to_path(f'{input_dir}').glob('*.vcf.gz')
        ]

    else:
        input_vcf_list = [input_dir]

    for vcf_file in input_vcf_list:
        vcf_input = b.read_input(vcf_file)
        if sharded:
            shard_num = (
                vcf_file.split('/')[-1].split('.')[0].split('_')[-1]
            )  # extracts shard number from file name
        else:
            shard_num = ''
        trtools_job = b.new_job(name=f'statSTR {caller} {shard_num}')
        trtools_job.image(TRTOOLS_IMAGE)
        trtools_job.storage(job_storage)
        trtools_job.memory(job_memory)

        trtools_job.declare_resource_group(ofile={'tab': '{root}.tab'})

        trtools_job.command(
            f"""
            set -ex;
            statSTR --vcf {vcf_input} --vcftype {caller} --out {trtools_job.ofile} --thresh --afreq --acount --hwep --het --entropy --mean --mode --var --numcalled

            """
        )
        if sharded:
            output_path_vcf = output_path(
                f'statSTR_shard_{shard_num}_{caller}', 'analysis'
            )
        else:
            output_path_vcf = output_path(f'statSTR_samples_{caller}', 'analysis')
        b.write_output(trtools_job.ofile, output_path_vcf)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
