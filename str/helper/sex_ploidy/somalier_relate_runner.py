#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script runs somalier relate

 analysis-runner --dataset "bioheart" \
    --description "Somalier relate runner" \
    --access-level "test" \
    --output-dir "hoptan-str/somalier" \
    somalier_relate_runner.py --input-dir=gs://cpg-bioheart-test/cram
t
"""

import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

SOMALIER_IMAGE = config['images']['somalier']


@click.option(
    '--input-dir',
    help='Input directory to cram.somalier files',
)
@click.option(
    '--job-storage', help='Storage of the Hail batch job eg 30G', default='10G'
)
@click.option('--job-memory', help='Memory of the Hail batch job', default='8G')
@click.option('--job-ncpu', help='Number of CPUs of the Hail batch job', default=4)
@click.command()
def main(
    input_dir,
    job_memory,
    job_ncpu,
    job_storage,
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()
    input_files = list(to_path(input_dir).glob('*.cram.somalier'))
    input_files = [
        str(gs_path) for gs_path in input_files
    ]  # coverts into a string type
    batch_input_files =[]
    for each_file in input_files:
        batch_input_files.append(
                        b.read_input(each_file))

    somalier_job = b.new_job(name=f'Somalier relate')
    somalier_job.image(SOMALIER_IMAGE)
    somalier_job.storage(job_storage)
    somalier_job.memory(job_memory)
    somalier_job.cpu(job_ncpu)

    somalier_job.declare_resource_group(
        somalier_output={
            'samples.tsv': '{root}.somalier_output.samples.tsv',
            'pairs.tsv': '{root}.somalier_output.pairs.tsv',
            'groups.tsv': '{root}.somalier_output.groups.tsv',
            '.html': '{root}.somalier_output.html',
        }
    )

    somalier_job.command(
        f"""
                somalier relate  \\
                {" ".join(batch_input_files)} \\
                --infer \\
                -o somalier_output
                """
    )
    # ExpansionHunter output writing
    somalier_output_path = output_path(f'somalier', 'analysis')
    b.write_output(somalier_job.somalier_output, somalier_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
