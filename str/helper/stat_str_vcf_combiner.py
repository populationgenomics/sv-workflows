#!/usr/bin/env python3

"""
This script combines sharded statSTR outputs into one output file, and assumes the genotyper used was ExpansionHunter.

analysis-runner --access-level test --dataset tob-wgs --description  \
    'stat STR combiner' --output-dir 'hoptan-str/shard_workflow_test/stat_str_vcf_combiner/v1' \
    stat_str_vcf_combiner.py \
    --input-dir=gs://cpg-tob-wgs-test/hoptan-str/shard_workflow_test/sharded_stat_str/sharded_stat_str
"""
import click

from cpg_utils import to_path
from cpg_utils.hail_batch import output_path


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
    input_file_paths = [str(gspath) for gspath in to_path(input_dir).glob('*.tab')]

    # Initialize variables to store information
    chrom_line = ''

    temporary_gt_file = 'temporary_gt_file.txt'
    with open(temporary_gt_file, 'w', encoding='utf-8') as handle:
        # Process each input file
        for index, input_file in enumerate(input_file_paths):
            input_file = to_path(input_file)
            print(f'Parsing {input_file}')

            with open(input_file, 'rt', encoding='utf-8') as f:
                for line in f:
                    # Collect information from the header lines
                    if line.startswith('chrom'):
                        if index == 1:
                            chrom_line = line
                    elif line.startswith('chr'):
                        # Collect calls after #CHROM in a temp file
                        handle.write(line)

    print(f'Parsed {len(input_file_paths)} sharded VCFs')


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
