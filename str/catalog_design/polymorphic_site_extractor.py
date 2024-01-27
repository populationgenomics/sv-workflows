#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This Hail Query script outputs a CSV file containing REPIDs (RepeatIDs) of sites that are polymorphic in a given VCF.
Polymorphic sites are those that have at least two distinct alleles occuring in the VCF.
Note this script also removes chrY and chrM sites.


 analysis-runner --dataset "bioheart" \
    --description "polymorphic-site-extractor" \
    --access-level "test" \
    --output-dir "hoptan-str/catalog-design/" \
    polymorphic_site_extractor.py --file-path=gs://cpg-tob-wgs-test/hoptan-str/shard_workflow_test/merge_str_vcf_combiner/combined_eh.vcf

"""

import hail as hl
import click
import csv

from cpg_utils import to_path
from cpg_utils.hail_batch import output_path, init_batch


def polymorphic_site_extractor(file_path, gcs_path):
    init_batch()
    # read in VCF into mt format
    hl.import_vcf(file_path).write('5M_n200.mt', overwrite=True)
    mt = hl.read_matrix_table('5M_n200.mt')

    mt = hl.sample_qc(mt)
    mt = hl.variant_qc(mt)

    ## remove chrY and chrM
    filtered_mt = mt.filter_rows(~hl.str(mt.locus.contig).startswith("chrY"))
    filtered_mt = filtered_mt.filter_rows(
        ~hl.str(filtered_mt.locus.contig).startswith("chrM")
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

    # Write the combined information to the output file
    with to_path(gcs_path).open('w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(rep_id_list)


@click.option(
    '--file-path',
    help='GCS file path VCF',
    type=str,
)
@click.command()
def main(file_path):
    gcs_output_path = output_path(f'polymorphic_sites_rep_id.csv', 'analysis')
    polymorphic_site_extractor(file_path, gcs_output_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
