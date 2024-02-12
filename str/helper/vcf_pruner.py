#!/usr/bin/env python3

"""
This script prunes a VCF file, removing all variants that are not present in a provided sharded catalog.
The output is a pruned VCF file, sharded in the same way as the input catalog, output to a GCS bucket.

Option to provide multiple samples and create pruned sharded VCFs for each sample.

analysis-runner --access-level test --dataset tob-wgs --description \
    'VCF pruner' --output-dir 'str/5M_run_combined_vcfs_pruned/v2' \
    vcf_pruner.py \
    --json-file-dir=gs://cpg-tob-wgs-test/str/5M_3M_merge_experiment/3M/CPG308296 \
    --vcf-file-dir=gs://cpg-tob-wgs-test/str/5M_3M_merge_experiment CPG308288
"""
import json
import click

from cpg_utils.config import get_config
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()


def variant_id_collector(catalog_file):
    """Collects all variant IDs from a sharded catalog file and returns a set of unique variant IDs."""

    variant_ids = []

    with to_path(catalog_file).open() as f:
        for line in f:
            if not line.startswith('#'):
                row_info = line.split('\t')[7]
                var_id = (row_info.split(';')[4])[6:]
                variant_ids.append(var_id)

    return set(variant_ids)

def pruner(vcf_file_path, cpg_id, chunk_number, variant_ids):
    """Prunes a VCF file, removing all variants that are not present in the list of provided variant IDs.
    Writes output to GCS bucket.
    """
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
                var_id = {(row_info.split(';')[4])[6:]}
                if var_id & variant_ids:
                    gt_lines.append(line)

    # Sort ALT lines alphabetically and convert to a list
    sorted_alt_lines = sorted(alt_lines)

    # Write the combined information to the output file
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
    """Prunes a VCF file, removing all variants that are not present in a provided sharded VCF."""
    # Initializing Batch
    b = get_batch()

    # list of catalog files (multiple, if catalog is sharded)
    catalog_files = list(to_path(json_file_dir).glob('*.vcf'))
    catalog_files = [
        str(gs_path) for gs_path in catalog_files
    ]  # coverts into a string type
    for catalog_file in catalog_files:
        variant_id_collector_job = b.new_python_job(
            name=f'Variant ID collector job: {catalog_file}'
        )
        variant_id_collector_job.memory('16G')
        variant_id_collector_job.storage('20G')
        variant_ids = variant_id_collector_job.call(variant_id_collector, catalog_file)

        for cpg_id in cpg_ids:
            # make input_files GSPath elements into a string type object
            vcf_file_path = f'{vcf_file_dir}/{cpg_id}_combined.vcf'
            chunk_number = catalog_file.split('/')[-1].split('_')[1].split('.')[0]
            vcf_pruner_job = b.new_python_job(
                name=f'VCF Combiner job: {cpg_id} chunk {chunk_number}'
            )
            vcf_pruner_job.memory('16G')
            vcf_pruner_job.storage('20G')
            vcf_pruner_job.call(
                pruner, vcf_file_path, cpg_id, chunk_number, variant_ids
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
