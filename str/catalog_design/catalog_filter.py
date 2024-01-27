#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script filters an ExpansionHunter JSON catalog file to only include loci with IDs provided in a separate CSV (eg polymorphic sites only)

 analysis-runner --dataset "tob-wgs" \
    --description "catalog-filter" \
    --access-level "test" \
    --output-dir "hoptan-str/catalog-design/" \
    catalog_filter.py --catalog-path=gs://cpg-tob-wgs-test/hoptan-str/5M_run/combined_catalog.trf_at_least_9bp.with_adjacent_loci.annotated_and_filtered.json \
    --filter-id-path=gs://cpg-tob-wgs-test/hoptan-str/5M_run/polymorphic_rep_id_n200.csv

"""
import pandas as pd
import csv
from itertools import chain
import json
import click


from cpg_utils import to_path
from cpg_utils.hail_batch import output_path, init_batch


def catalog_filter(catalog_path, filter_id_path, gcs_path):
    polymorphic_rep_id =[]

    # Read the CSV file and append each row to the list
    with to_path(filter_id_path).open('r') as file:
        reader = csv.reader(file)
        polymorphic_rep_id = map(list, zip(*reader))
    print(f'Parsed {len(polymorphic_rep_id)} polymorphic REPIDs')

    flattened_polymorphic_rep_id = list(chain.from_iterable(polymorphic_rep_id))
    polymorphic_rep_id_set =set(flattened_polymorphic_rep_id)

    with to_path(catalog_path).open('r') as json_file:
        # Load the JSON content
        catalog = json.load(json_file)

    # Filter the JSON entries based on the VariantID
    filtered_data =[]

    for entry in catalog:
        # Check if 'VariantId' exists, use 'LocusId' otherwise
        entry_variant_ids = set(entry['VariantId']) if 'VariantId' in entry else {entry['LocusId']}

        # Check if there is an intersection with polymorphic_rep_id_set
        if entry_variant_ids & polymorphic_rep_id_set:
            filtered_data.append(entry)

    # Write the combined information to the output file
    with to_path(gcs_path).open('w') as out_file:
        json.dump(filtered_data, out_file, indent=2)

@click.option(
    '--filter-id-path',
    help='GCS file path VCF',
    type=str,
)
@click.option(
    '--catalog-path',
    help='GCS file path to catalog JSON',
    type=str,
)
@click.command()
def main(filter_id_path, catalog_path):
    gcs_output_path = output_path(f'filtered_catalog_json', 'analysis')
    catalog_filter(catalog_path, filter_id_path, gcs_output_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter