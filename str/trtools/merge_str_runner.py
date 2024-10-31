#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals
"""
This script merges ExpansionHunter vcf.gz files into one combined VCF.
Please ensure merge_prep.py has been run on the vcf files prior to running mergeSTR.py

Optional ability to add in VCFs from another file directory (but must be sharded in the same way)
Specify VCFs from each distinct file directory as a comma separated list of the input directory and the sample list file (in that order) eg: input-dir-1,sample-1 input-dir-2,sample-2

For example:
analysis-runner --access-level full --dataset tob-wgs --description '5M-3M mergeSTR tester' --output-dir 'str/5M-3M experiment/merge_str/v1' merge_str_runner.py --num-shards=27 \
gs://cpg-tob-wgs-main-analysis/str/5M_run_combined_vcfs_pruned/merge_str_prep/v4,gs://cpg-tob-wgs-test/str/polymorphic_run/mergeSTR-tester-5M.txt \
gs://cpg-tob-wgs-main-analysis/str/polymorphic_run/merge_str_prep/v1,gs://cpg-tob-wgs-test/str/polymorphic_run/mergeSTR-tester-3M.txt

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os

import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

TRTOOLS_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/trtools:v6.0.2'


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
@click.option('--job-cpu', help='Number of CPUs of the Hail batch job', default=1)
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
        cpg_ids = []
        for pair in input_file_paths:
            input_dir, sample_list = pair.split(',')
            with to_path(sample_list).open() as f:
                ids = [line.strip() for line in f]
                cpg_ids.extend(ids)
                for id in ids:
                    each_vcf = os.path.join(input_dir, f'{id}_eh_shard{shard_index}.reheader.vcf.gz')
                    batch_vcfs.append(
                        '${BATCH_TMPDIR}/'+str(b.read_input_group(
                            **{
                                'vcf.gz': each_vcf,
                                'vcf.gz.tbi': f'{each_vcf}.tbi',
                            },
                        )['vcf.gz'].split('/', 1)[-1])
                    )
                num_samples = num_samples + len(ids)

        if len(cpg_ids) != len(set(cpg_ids)):
            raise ValueError('Duplicate CPG IDs detected in sample list')

        batch_vcfs_list = '${BATCH_TMPDIR}/batch_vcfs.list'

        newline = '\n'
        # Writing cat...EOF on one line avoids any indentation on the batch_vcfs and EOF lines
        trtools_job.command(
            f"""
        cat <<'EOF' >{batch_vcfs_list}\n{newline.join(batch_vcfs)}\nEOF

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
