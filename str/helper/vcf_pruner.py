#!/usr/bin/env python3

"""
This script prunes a VCF file (vcf-file-dir), removing all variants that are not present in a separate provided sharded VCF (vcf-catalog-dir).
This script assumes that the input VCF file to be pruned contains variants that are a super-set of the variants in vcf-catalog-dir.

The output is a pruned VCF file, sharded in the same way as the files in vcf-catalog-dir, output to a GCS bucket.

Optional to provide multiple samples and create pruned sharded VCFs for each sample.

analysis-runner --access-level test --dataset tob-wgs --description \
    'VCF pruner' --output-dir 'str/5M_run_combined_vcfs_pruned/v2' \
    vcf_pruner.py \
    --vcf-catalog-dir=gs://cpg-tob-wgs-test/str/5M_3M_merge_experiment/3M/CPGXXX \
    --vcf-file-dir=gs://cpg-tob-wgs-test/str/5M_3M_merge_experiment CPGXXXXX
"""
import click

from cpg_utils.config import get_config
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()


def variant_id_collector(catalog_file):
    """Collects all variant IDs from a sharded VCF file and returns a list of variant IDs."""

    variant_ids = []

    with to_path(catalog_file).open() as f:
        for line in f:
            if not line.startswith('#'):
                # wrangle variant ID from INFO column, assumes ExpansionHunter VCF format
                row_info = line.split('\t')[7]
                var_id = (row_info.split(';')[4])[6:]

                variant_ids.append(var_id)

    return variant_ids


def pruner(vcf_file_path, cpg_id, chunk_number, variant_id_order):
    """Prunes a VCF file, removing all variants that are not present in the list of provided variant IDs.
    Writes output to GCS bucket.
    """
    # Initialize variables to store information
    fileformat_line = ''
    info_lines = []
    alt_lines = set()
    chrom_line = ''
    gt_lines = {}  # key = variant ID, value = line in VCF containing GT
    # set structure allows for faster intersection
    variant_id_set = set(variant_id_order)

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

                if (
                    var_id & variant_id_set
                ):  # if the variant ID is in the set of target variant IDs
                    var_id = ''.join(
                        map(str, var_id)
                    )  # converts var_id from set to string
                    gt_lines[var_id] = line

    # Sort ALT lines alphabetically and convert to a list
    sorted_alt_lines = sorted(alt_lines)

    # sort gt_lines by variant ID in the order provided by variant_id_order
    sorted_gt = {key: gt_lines[key] for key in variant_id_order if key in gt_lines}
    print(f'Parsed {len(sorted_gt)} variants from {cpg_id} shard {chunk_number}')

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
        out_file.writelines(sorted_gt.values())


@click.command()
@click.option(
    '--vcf-catalog-dir',
    help='Parent input directory for sharded VCFs containing target variants (subfolders should be labelled with CPG ID)',
)
@click.option(
    '--vcf-file-dir',
    help='GCS path to folder containg the VCF file(s) to be pruned',
)
@click.argument('cpg-ids', nargs=-1)
def main(vcf_catalog_dir, vcf_file_dir, cpg_ids: list[str]):
    """Prunes a VCF file, removing all variants that are not present in a provided sharded VCF."""
    # Initializing Batch
    b = get_batch()

    # list of catalog files (multiple, if catalog is sharded)
    catalog_files = list(to_path(vcf_catalog_dir).glob('*.vcf'))
    catalog_files = [
        str(gs_path) for gs_path in catalog_files
    ]  # coverts into a string type
    for catalog_file in catalog_files:
        variant_id_collector_job = b.new_python_job(
            name=f'Variant ID collector job: {catalog_file}'
        )
        variant_id_collector_job.memory('8G')
        variant_id_collector_job.storage('10G')
        variant_id_order = variant_id_collector_job.call(
            variant_id_collector, catalog_file
        )

        for cpg_id in cpg_ids:

            vcf_file_path = f'{vcf_file_dir}/{cpg_id}_combined.vcf'
            # extract shard number from the VCF file name
            chunk_number = catalog_file.split('/')[-1].split('shard')[1].split('.')[0]
            vcf_pruner_job = b.new_python_job(
                name=f'VCF Combiner job: {cpg_id} chunk {chunk_number}'
            )
            vcf_pruner_job.memory('8G')
            vcf_pruner_job.storage('10G')
            vcf_pruner_job.call(
                pruner,
                vcf_file_path,
                cpg_id,
                chunk_number,
                variant_id_order,
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
