#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script:
(1) parses a given VCF for polymorphic sites
(2) filters a given catalog JSON file to only include those sites, outputting the filtered catalog to a new JSON file.
(3) shards the filtered catalog into chunks, writing out each chunk as a separate JSON file.

Polymorphic sites are those that have at least two distinct alleles occuring in the VCF.
Note this script filters out chrY and chrM sites by default.

 analysis-runner --dataset "tob-wgs" \
    --description "polymorphic-site-extractor" \
    --access-level "test" \
    --output-dir "hoptan-str/catalog-design/" \
    --memory=32G --storage=20G \
    polymorphic_site_extractor.py --vcf-path=gs://cpg-tob-wgs-test/hoptan-str/shard_workflow_test/merge_str_vcf_combiner/combined_eh.vcf \
    --catalog-path=gs://cpg-tob-wgs-test/hoptan-str/5M_run/5M_sharded_100k/chunk_1.json
    --folder_name=sharded_polymorphic_catalog

"""
import json

import click

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import init_batch, output_path


def polymorphic_site_extractor(file_path):
    """Extracts polymorphic sites from a VCF file and returns a list of REPIDs (similar to rsids) representing those sites
    The output has been benchmarked with Gymrek's statSTR and agrees with the number of polymorphic sites found
    """
    init_batch()
    # read in VCF into mt format
    mt = hl.import_vcf(file_path)

    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)

    # remove chrY and chrM, monomorphic REF sites, monomorphic ALT at bialellic loci
    filtered_mt = mt.filter_rows(
        (hl.str(mt.locus.contig).startswith('chrY'))
        | (hl.str(mt.locus.contig).startswith('chrM'))
        | (hl.len(mt.alleles) == 1)
        | ((hl.len(mt.variant_qc.AC) == 2) & (mt.variant_qc.AC[0] == 0)),
        keep=False,
    )

    # collect the REPIDs into one list
    rep_id_list = filtered_mt.info.REPID.collect()
    rep_id_set = set(rep_id_list)

    return rep_id_set


def catalog_filter(polymorphic_rep_id_set: set[str], catalog_path: str, gcs_output_path: str):
    """Retains loci in a JSON file that intersect with a list of REPIDs (STR equiv. of rsids) and writes the filtered catalog to a new JSON file"""
    with to_path(catalog_path).open('r') as json_file:
        # Load the JSON content
        catalog = json.load(json_file)

    filtered_data = []

    for entry in catalog:
        # Check if 'VariantId' exists, use 'LocusId' otherwise
        entry_variant_ids = set(entry['VariantId']) if 'VariantId' in entry else {entry['LocusId']}

        # Append to filtered_data if the entry Variant Id intersects with polymorphic_rep_id_set
        if entry_variant_ids & polymorphic_rep_id_set:
            filtered_data.append(entry)

    # Write to output
    with to_path(gcs_output_path).open('w') as out_file:
        json.dump(filtered_data, out_file, indent=2)
    return filtered_data


def catalog_sharder(filtered_data, chunk_size, folder_name):
    """Shards a filtered catalog JSON file into chunks of size chunk_size"""
    if not isinstance(filtered_data, list):
        raise ValueError('Invalid JSON format. The file should contain a list.')

    total_entries = len(filtered_data)
    num_chunks = total_entries // chunk_size
    remainder = total_entries % chunk_size

    if remainder > 0:
        num_chunks += 1

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = (i + 1) * chunk_size
        chunk_data = filtered_data[start_idx:end_idx]

        output_file_path = output_path(f'{folder_name}/chunk_{i + 1}.json', 'analysis')

        with to_path(output_file_path).open('w') as output_file:
            json.dump(chunk_data, output_file, indent=2)

        print(f'Chunk {i + 1} created: {output_file_path}')

    print(f'Total entries: {total_entries}')
    print(f'{num_chunks} chunks written successfully!')


@click.option('--vcf-path', help='GCS file path VCF')
@click.option('--catalog-path', help='GCS file path to catalog JSON to be filtered')
@click.option(
    '--chunk-size',
    help='Number of entries per shard',
    type=int,
    default=100000,
)
@click.option(
    '--folder-name',
    help='Name of output folder to store shards',
    default='sharded_polymorphic_catalog',
)
@click.command()
def main(vcf_path, catalog_path, chunk_size, folder_name):
    # Extract polymorphic sites from VCF
    rep_id_set = polymorphic_site_extractor(vcf_path)

    # Filter catalog JSON file for polymorphic sites
    catalog_output_path = output_path('filtered_polymorphic_catalog.json', 'analysis')
    filtered_catalog = catalog_filter(rep_id_set, catalog_path, catalog_output_path)

    # Shard the filtered catalog
    catalog_sharder(filtered_catalog, chunk_size, folder_name)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
