#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This Hail Query script parses a given VCF for polymorphic sites and then filters a given catalog JSON file to only include those sites.
Polymorphic sites are those that have at least two distinct alleles occuring in the VCF.
Note this script also removes chrY and chrM sites.

 analysis-runner --dataset "tob-wgs" \
    --description "polymorphic-site-extractor" \
    --access-level "test" \
    --output-dir "hoptan-str/catalog-design/" \
    --memory=32G --storage=20G \
    polymorphic_site_extractor.py --vcf-path=gs://cpg-tob-wgs-test/hoptan-str/shard_workflow_test/merge_str_vcf_combiner/combined_eh.vcf \
    --catalog-path=gs://cpg-tob-wgs-test/hoptan-str/5M_run/combined_catalog.trf_at_least_9bp.with_adjacent_loci.annotated_and_filtered.json

"""
import json
import hail as hl
import click


from cpg_utils import to_path
from cpg_utils.hail_batch import output_path, init_batch


def polymorphic_site_extractor(file_path):
    """ Extracts polymorphic sites from a VCF file and returns a list of REPIDs (similar to rsids) representing those sites """
    init_batch()
    # read in VCF into mt format
    mt = hl.import_vcf(file_path)

    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)

    # remove chrY and chrM
    filtered_mt = mt.filter_rows(~hl.str(mt.locus.contig).startswith('chrY'))
    filtered_mt = filtered_mt.filter_rows(
        ~hl.str(filtered_mt.locus.contig).startswith('chrM')
    )

    # remove loci that are monomorphic for the REF allele
    filtered_mt = filtered_mt.filter_rows(hl.len(filtered_mt.alleles) > 1)

    # at a biallelic locus, remove loci that are monomorphic for 1 ALT allele
    filtered_mt = filtered_mt.filter_rows(
        ~(
            (hl.len(filtered_mt.variant_qc.AC) == 2)
            & (filtered_mt.variant_qc.AC[0] == 0)
        )
    )

    # collect the REPIDs into one list
    rep_id_list = filtered_mt.info.REPID.collect()

    return rep_id_list


def catalog_filter(rep_id_list, catalog_path, gcs_output_path):
    """ Retains loci in a JSON file that intersect with a list of REPIDs (rsids) and writes the filtered catalog to a new JSON file"""
    with to_path(catalog_path).open('r') as json_file:
        # Load the JSON content
        catalog = json.load(json_file)

    # Filter the JSON entries based on the VariantID
    filtered_data = []

    # Make the rep_id_list into a set
    polymorphic_rep_id_set = set(rep_id_list)

    for entry in catalog:
        # Check if 'VariantId' exists, use 'LocusId' otherwise
        entry_variant_ids = (
            set(entry['VariantId']) if 'VariantId' in entry else {entry['LocusId']}
        )

        # Check if there is an intersection with polymorphic_rep_id_set
        if entry_variant_ids & polymorphic_rep_id_set:
            filtered_data.append(entry)

    # Write the combined information to the output file
    with to_path(gcs_output_path).open('w') as out_file:
        json.dump(filtered_data, out_file, indent=2)


@click.option(
    '--vcf-path',
    help='GCS file path VCF',
    type=str,
)
@click.option(
    '--catalog-path',
    help='GCS file path to catalog JSON',
    type=str,
)
@click.command()
def main(vcf_path, catalog_path):
    gcs_output_path = output_path(f'filtered_polymorphic_catalog.json', 'analysis')
    rep_id_list = polymorphic_site_extractor(vcf_path)
    catalog_filter(rep_id_list, catalog_path, gcs_output_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
