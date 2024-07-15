#!/usr/bin/env python3

"""
This script prepares two files:

- Matrix of beta and SE for eSTRs passing FDR <5% across all cell types

- Matrix of randomly sampled beta and SE for 20k eSTRs.

analysis-runner --dataset "bioheart" \
    --description "Prepare inputs for mashr" \
    --access-level "test" \
    --output-dir "potato" \
    prep_inputs.py
"""

import os
import sys

import numpy as np
import pandas as pd

from cpg_utils.hail_batch import get_batch


def cell_chrom_parser(cell, chrom, estrs_coord_chrom):
    cell_df = pd.DataFrame()
    for index, row in estrs_coord_chrom.iterrows():
        chrom = row['chr']
        gene_name = row['gene_name']
        pos = row['pos']
        motif = row['motif']
        ref_len = row['ref_len']
        try:
            df = pd.read_csv(
                f'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/{cell}/{chrom}/{gene_name}_100000bp_meta_results.tsv',
                sep='\t',
            )
            df = df[(df['pos'] == pos) & (df['motif'] == motif) & (df['ref_len'] == ref_len)]
        except:
            continue

        beta = df['coeff_meta'].iloc[0]
        se = df['se_meta'].iloc[0]
        # save the beta and se into a row:
        beta_se = pd.DataFrame(
            {
                'chrom': chrom,
                'pos': pos,
                'motif': motif,
                'ref_len': ref_len,
                'gene': gene_name,
                f'{cell}_beta': beta,
                f'{cell}_se': se,
            },
            index=[0],
        )
        cell_df = pd.concat([cell_df, beta_se], axis=0)
    cell_df.to_csv(
        'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/estrs_beta_se/{cell}/{chrom}/beta_se.tsv',
        sep='\t',
        index=False,
    )


def main():
    b = get_batch(name='Prep eSTRs for mashr')
    cell_types = 'CD4_TCM'
    celltypes = cell_types.split(',')
    # load in the list of eSTRs passing FDR <5% across all cell types:
    estrs_coord = pd.read_csv(
        'gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/mashr/estrs_coord_gene.csv',
    )
    for cell in celltypes:
        for chrom in range(1, 2):
            estrs_coord_chrom = estrs_coord[estrs_coord['chr'] == f'chr{chrom}']
            job = b.new_python_job(f'Prep eSTRs for mashr {cell} {chrom}')
            job.call(cell_chrom_parser, cell, chrom, estrs_coord_chrom)
    b.run(wait=False)


if __name__ == '__main__':
    main()
