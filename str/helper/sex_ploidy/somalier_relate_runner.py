#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script runs somalier relate on a set of cram.somalier files provided as input directory file path.

 analysis-runner --dataset "bioheart" \
    --description "Somalier relate runner" \
    --access-level "standard" \
    --output-dir "qc-stand-alone/tob_bioheart/somalier" \
    somalier_relate_runner.py --input-dir-1=gs://cpg-bioheart-main/cram --input-dir-2=gs://cpg-tob-wgs-main/cram

"""

import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

SOMALIER_IMAGE = config['images']['somalier']


@click.option(
    '--input-dir-1',
    help='Input directory to cram.somalier files',
)
@click.option(
    '--input-dir-2',
    help='2nd (optional)input directory to cram.somalier files',
    default=None,
)
@click.option(
    '--job-storage', help='Storage of the Hail batch job eg 30G', default='10G'
)
@click.option('--job-memory', help='Memory of the Hail batch job', default='8G')
@click.option('--job-ncpu', help='Number of CPUs of the Hail batch job', default=4)
@click.command()
def main(
    job_memory, job_ncpu, job_storage, input_dir_1, input_dir_2
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()
    input_files = list(to_path(input_dir_1).glob('*.cram.somalier'))
    input_files = [
        str(gs_path) for gs_path in input_files
    ]  # coverts into a string type
    if input_dir_2 is not None:
        input_files_2 = list(to_path(input_dir_2).glob('*.cram.somalier'))
        input_files_2 = [str(gs_path) for gs_path in input_files_2]
        input_files.extend(input_files_2)

    num_samples = len(input_files)
    batch_input_files = []
    for each_file in input_files:
        batch_input_files.append(b.read_input(each_file))

    somalier_job = b.new_job(name=f'Somalier relate: {num_samples} samples')
    somalier_job.image(SOMALIER_IMAGE)
    somalier_job.storage(job_storage)
    somalier_job.memory(job_memory)
    somalier_job.cpu(job_ncpu)

    somalier_job.command(
        f"""
                somalier relate  \\
                {" ".join(batch_input_files)} \\
                --infer \\
                -o related

                mv related.pairs.tsv {somalier_job.output_pairs}
                mv related.samples.tsv {somalier_job.output_samples}
                mv related.html {somalier_job.output_html}
                """
    )
    # Output writing
    b.write_output(
        somalier_job.output_samples,
        str(output_path(f'{num_samples}_samples_somalier.samples.tsv', 'analysis')),
    )

    b.write_output(
        somalier_job.output_pairs,
        str(output_path(f'{num_samples}_samples_somalier.pairs.tsv', 'analysis')),
    )
    b.write_output(
        somalier_job.output_html,
        str(output_path(f'{num_samples}_samples_somalier.html', 'analysis')),
    )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
