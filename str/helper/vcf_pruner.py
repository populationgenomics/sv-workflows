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
    '--json-file-dir',
    help='Parent input directory for sharded VCFs (subfolders should be labelled with CPG ID)',
)
@click.option(
    '--vcf-file-dir',
    help='GCS path to folder containg the VCF file(s) to be pruned',
)
@click.argument('cpg-ids', nargs=-1)

def main(json_file_dir, vcf_file_dir, cpg_ids: list[str]):

    # list of catalog files (multiple, if catalog is sharded)
    catalog_files = list(to_path(json_file_dir).glob('*.json'))
    catalog_files = [
        str(gs_path) for gs_path in catalog_files
    ]  # coverts into a string type
    for catalog_file in catalog_files:
        variant_ids= []

        with to_path(catalog_file).open('r') as json_file:
            # Load the JSON content
            catalog = json.load(json_file)

            for entry in catalog:
                # Check if 'VariantId' exists, use 'LocusId' otherwise
                entry_variant_ids = (
                    entry['VariantId'] if 'VariantId' in entry else [entry['LocusId']]
                )
                variant_ids.extend(entry_variant_ids)
        variant_ids = set(variant_ids)


        # Initialize variables to store information
        fileformat_line = ''
        info_lines = []
        alt_lines = set()
        chrom_line = ''
        gt_lines = []

        for cpg_id in cpg_ids:
            # make input_files GSPath elements into a string type object
            vcf_file_path = f'{vcf_file_dir}/{cpg_id}_combined.vcf'
            # Process each input file

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
                        var_id = {(row_info.split(';')[4])[6:]}
                        if var_id & variant_ids:
                            gt_lines.append(line)

            # Sort ALT lines alphabetically and convert to a list
            sorted_alt_lines = sorted(alt_lines)

            # Write the combined information to the output file
            chunk_number = catalog_file.split('/')[-1].split('_')[1].split('.')[0]
            gcs_out_path = output_path(f'{cpg_id}/{cpg_id}_eh_shard{chunk_number}.vcf')
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
