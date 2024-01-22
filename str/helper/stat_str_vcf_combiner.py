#!/usr/bin/env python3

"""
This script combines sharded statSTR outputs into one output file, and assumes the genotyper used was ExpansionHunter.

analysis-runner --access-level test --dataset tob-wgs --description  \
    'stat STR combiner' --output-dir 'hoptan-str/shard_workflow_test/stat_str_vcf_combiner/v1' \
    stat_str_vcf_combiner.py \
    --input-dir=gs://cpg-tob-wgs-test/hoptan-str/shard_workflow_test/sharded_stat_str
"""
import gzip
import click

from cpg_utils import to_path
from cpg_utils.hail_batch import output_path

# custom sorted order of shards - PLEASE UPDATE FOR LATER USE
SORTED_KEY_ORDER = [
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
    # 46,
]


@click.command()
@click.option('--input-dir', help='Parent input directory for sharded statSTR outputs')
@click.option('--output', help='Name of output VCF', default='statSTR_combined.tab')
def main(input_dir, output):
    """
    Combines sharded statSTR output tabs in input_dir into one combined statSTR tab file,
    writing it to a GCS output path
    Takes an input directory containing statSTR sharded outputs
    Aggregates all sharded data into a single output file

    Doesn't use batch any more, as mem, cpu and disk can all be analysis-runner CLI args
    """

    # make input_files GSPath elements into a string type object
    input_file_paths = map(str, to_path(input_dir).glob('*.tab'))

    # create a dictionary where key is the shard number and value is the path to the shard
    input_files_dict = {
        int(file_path.split('_shard_shard')[1].split('_eh')[0]): file_path
        for file_path in input_file_paths
    }
    print(input_files_dict.keys())

    # Initialize variables to store information
    chrom_line = ''

    temporary_gt_file = 'temporary_gt_file.txt'
    with open(temporary_gt_file, 'w', encoding='utf-8') as handle:
        # Process each input file
        for key in SORTED_KEY_ORDER:
            input_file = to_path(input_files_dict[key])
            print(f'Parsing {input_file}')

            with gzip.open(input_file, 'rt') as f:
                for line in f:
                    # Collect information from the header lines
                    if line.startswith('chrom'):
                        if key == 1:
                            chrom_line = line
                    elif line.startswith('chr'):
                        # Collect calls after #CHROM in a temp file
                        handle.write(line)

    print(f'Parsed {len(SORTED_KEY_ORDER)} sharded VCFs')

    # Write the combined information to the output file
    with to_path(output_path(output, 'analysis')).open('w') as out_file:
        # Write CHROM line
        out_file.write(chrom_line)

        # read-write all GT lines from temporary file
        with open(temporary_gt_file, 'r', encoding='utf-8') as handle:
            for line in handle:
                out_file.write(line)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
