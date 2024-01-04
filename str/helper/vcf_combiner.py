#!/usr/bin/env python3
# pylint: disable=import-error, too-many-locals

"""
This script combines sharded VCFs (from one sample and one caller) into one VCF.

analysis-runner --access-level test --dataset tob-wgs --description \
    'VCF combiner' --output-dir 'str/5M_run_combined_vcfs' \
    vcf_combiner.py \
    --input-dir=gs://cpg-tob-wgs-test/hoptan-str/sharded_tester
"""
import click

from cpg_utils.config import get_config
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()


def combine_vcf_files(input_dir, gcs_out_path):
    """Combines sharded VCFs in input_dir into one combined VCF, writing it to a GCS output path"""
    input_files = to_path(input_dir).glob('*.vcf')

    # Initialize variables to store information
    fileformat_line = ''
    info_lines = []
    alt_lines = set()
    chrom_line = ''
    gt_lines = []

    # Process each input file
    for file_index, input_file in enumerate(input_files):
        with to_path(input_file).open() as f:
            for line in f:
                # Collect information from the header lines
                if line.startswith('##fileformat') and file_index == 0:
                    fileformat_line = line
                elif (
                    line.startswith('##INFO')
                    or line.startswith('##FILTER')
                    or line.startswith('##FORMAT')
                ) and file_index == 0:
                    info_lines.append(line)
                elif line.startswith('##ALT'):
                    # Collect ALT lines from all files into a set to remove duplicates
                    alt_lines.add(line)
                elif line.startswith('#CHROM') and file_index == 0:
                    chrom_line = line
                else:
                    # Collect calls after #CHROM
                    gt_lines.append(line)

    # Sort ALT lines alphabetically and convert to a list
    sorted_alt_lines = sorted(alt_lines)

    # Write the combined information to the output file
    with to_path(gcs_out_path).open('w') as out_file:
        # Write fileformat line
        out_file.write(fileformat_line)
        # Write INFO, FILTER, and FORMAT lines
        out_file.writelines(info_lines)
        # Write ALT lines
        out_file.writelines(sorted_alt_lines)
        # Write CHROM line
        out_file.write(chrom_line)
        # Write GT or lines containing the calls
        out_file.writelines(gt_lines)


@click.command()
@click.option('--input-dir', help='Input directory for sharded VCFs')
def main(input_dir):
    # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()
    cpg_id = input_dir.split('/')[-1]

    combiner_job = b.new_python_job(name=f'VCF Combiner job: {cpg_id}')

    out_path = output_path(f'{cpg_id}_combined.vcf', 'analysis')
    combiner_job.call(combine_vcf_files, input_dir, out_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
