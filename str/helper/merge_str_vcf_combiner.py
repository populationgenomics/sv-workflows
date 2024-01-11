#!/usr/bin/env python3

"""
This script combines sharded output VCFs from mergeSTR into one VCF, and assumes the genotyper used was ExpansionHunter.

analysis-runner --access-level test --dataset tob-wgs --description \
    'VCF combiner' --output-dir 'hoptan-str/shard_workflow_test/merge_str_vcf_combiner' \
    merge_str_vcf_combiner.py \
    --input-dir=gs://cpg-tob-wgs-test/hoptan-str/shard_workflow_test/merge_str
"""
import gzip
import click

from cpg_utils.config import get_config
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()


def combine_vcf_files(input_dir, gcs_out_path):
    """Combines sharded mergeSTR output VCFs in input_dir into one combined VCF, writing it to a GCS output path"""

    input_file_paths = to_path(input_dir).glob('*.vcf.gz')

    # make input_files GSPath elements into a string type object
    input_file_paths = {str(gs_path) for gs_path in input_file_paths}

    # create a dictionary where key is the shard number and value is the path to the shard
    input_files_dict = {}
    for file_path in input_file_paths:
        key = int(file_path.split('eh_shard')[1].split('.')[0])
        input_files_dict[key] = file_path

    # custom sorted order of shards - PLEASE UPDATE FOR LATER USE
    sorted_key_order = [
        1,
        12,
        23,
        34,
        45,
        47,
        48,
        49,
        50,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        32,
        33,
        35,
        36,
        37,
        38,
        39,
        40,
        41,
        42,
        43,
        44,
        46,
    ]

    # Initialize variables to store information
    fileformat_line = ''
    info_lines = []
    alt_lines = set()
    chrom_line = ''
    gt_lines = []

    shard_counter = 0
    # Process each input file
    for shard_counter, key in enumerate(sorted_key_order):
        input_file = to_path(input_files_dict[key])
        print(f'Parsing {input_file}')

        with gzip.open(input_file, 'rt') as f:
            for line in f:
                # Collect information from the header lines
                if line.startswith('##fileformat'):
                    if key == 1:  # first file processed is shard_1
                        fileformat_line = line
                elif (
                    line.startswith('##INFO')
                    or line.startswith('##FILTER')
                    or line.startswith('##FORMAT')
                    or line.startswith('##contig')
                    or line.startswith('##command')
                ):
                    if key == 1:
                        info_lines.append(line)
                elif line.startswith('##ALT'):
                    # Collect ALT lines from all files into a set to remove duplicates
                    alt_lines.add(line)
                elif line.startswith('#CHROM'):
                    if key == 1:
                        chrom_line = line
                elif not line.startswith('#'):
                    # Collect calls after #CHROM
                    gt_lines.append(line)

    print(f'Parsed {len(sorted_key_order)} sharded VCFs')
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
@click.option(
    '--input-dir',
    help='Parent input directory for sharded VCFs',
)
def main(input_dir):
    """
    Takes an input directory containing vcf shards
    Aggregates all sharded data into a single output file
    """

    out_path = output_path(f'combined_eh.vcf', 'analysis')
    combine_vcf_files(input_dir, out_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
