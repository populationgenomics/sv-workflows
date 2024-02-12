#!/usr/bin/env python3

"""
This script prunes a VCF file, removing all variants that are not present in a provided sharded catalog.
The output is a pruned VCF file, sharded in the same way as the input catalog, output to a GCS bucket.

analysis-runner --access-level test --dataset tob-wgs --description \
    'VCF combiner' --output-dir 'str/5M_run_combined_vcfs' \
    vcf_combiner.py \
    --input-dir=gs://cpg-tob-wgs-test/hoptan-str \
    sharded_tester
"""

import click
import json

from cpg_utils.config import get_config
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

@click.command()
@click.option(
    '--json-file-path',
    help='Parent input directory for sharded VCFs (subfolders should be labelled with CPG ID)',
)
@click.option(
    '--vcf-file-path',
    help='GCS path to the VCF file to be pruned',
)

def main(json_file_path, vcf_file_path):
    with to_path(json_file_path).open('r') as json_file:
        # Load the JSON content
        catalog = json.load(json_file)

    for entry in catalog:
        # Check if 'VariantId' exists, use 'LocusId' otherwise
        entry_variant_ids = (
            set(entry['VariantId']) if 'VariantId' in entry else {entry['LocusId']}
        )
    # Initialize variables to store information
    fileformat_line = ''
    info_lines = []
    alt_lines = set()
    chrom_line = ''
    gt_lines = []

    with to_path(vcf_file_path).open() as f:
        for line in f:
            # Collect information from the header lines
            if line.startswith('##fileformat'):
                fileformat_line = line
            elif (
                line.startswith('##INFO')
                or line.startswith('##FILTER')
                or line.startswith('##FORMAT')
            ):

                info_lines.append(line)
            elif line.startswith('##ALT'):
                # Collect ALT lines from all files into a set to remove duplicates
                alt_lines.add(line)
            elif line.startswith('#CHROM'):
                chrom_line = line
            elif not line.startswith('#'):
                # Collect calls after #CHROM
                row_info = line.split('\t')[7]
                var_id = row_info.split(';')[4][6:]
                if var_id in entry_variant_ids:
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
