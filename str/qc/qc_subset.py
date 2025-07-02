#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script annotates the ExpansionHunter MT with the following annotations:
- rep_length_alleles: repeat length of each allele
- motif_length: length of the motif
- bp_length_alleles: bp length of each allele
- allele_1_rep_length: repeat length of allele 1
- allele_2_rep_length: repeat length of allele 2
- allele_1_bp_length: bp length of allele 1
- allele_2_bp_length: bp length of allele 2
- aggregated_info.mode_allele:mode allele at each locus
- num_alleles: number of alleles at each locus
- allele_1_minus_mode: difference between allele 1 and mode allele
- allele_2_minus_mode: difference between allele 2 and mode allele
- sum_alleles_is_not_mode: count of non-mode alleles at each locus
- prop_alleles_is_not_mode: proportion of alleles that are not the mode allele per locus
- binom_hwep: binomial Hardy-Weinberg equilibrium p-value
- obs_het: proportion of observed heterozygous calls per locus

analysis-runner --access-level "test" --dataset "bioheart" --description "QC annotator" --output-dir "str/polymorphic_run/mt/bioheart_tob/n975_bioheart" qc_subset.py \
--mt-path=gs://cpg-bioheart-test/str/polymorphic_run/mt/bioheart_tob/v1_n1925/str_annotated.mt

"""

import click

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch, output_path
import pandas as pd

config = get_config()


@click.option('--mt-path', help='GCS Path to the input MT')
@click.command()
def main(mt_path):
    """
    Annotates the ExpansionHunter MT, and outputs annotated MT to GCS
    """

    init_batch(worker_memory='highmem')
    mt = hl.read_matrix_table(mt_path)
    print(f'MT dimensions: {mt.count()}')

    bioheart_ids = pd.read_csv('gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/bioheart_n975_sample_covariates.csv')['sample_id']
    tob_ids = pd.read_csv('gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/tob_n950/covariates/6_rna_pcs/CD4_TCM_covariates.csv')['sample_id']
    samples = bioheart_ids.to_list()

    # filter the MT to only include samples in the sample list
    mt = mt.filter_cols(hl.literal(samples).contains(mt.s))


    # write out mt to GCS path
    mt.write(output_path('str_annotated.mt'))

    # print mt schema
    mt.describe()


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
