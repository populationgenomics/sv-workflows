#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

import json

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

def main():
    dfs=[]
    genes = pd.read_csv('gs://cpg-bioheart-test/str/CD4_TCM_lead_TR_genes.csv')
    for row in genes.iterrows():
        gene = row['gene']
        chromosome = row['chrom']
        eqtl_results = pd.read_csv(
                f'gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results/CD4_TCM/{chromosome}/{gene}_100000bp_meta_results.tsv',
                sep='\t',
            )
        # get row(s) with minimum p-value
        min_pval = eqtl_results['pval_meta'].min()
        smallest_pval_rows = eqtl_results[eqtl_results['pval_meta'] == min_pval]
        #get coordinates of the lead TR
        lead_tr_coord = str(chromosome) + ':' + str(smallest_pval_rows.iloc[0]['pos'])
        lead_tr_motif = smallest_pval_rows.iloc[0]['motif']

        results_df = pd.DataFrame(
                {

                    'gene': [gene],
                    'lead_tr_coord': [lead_tr_coord],
                    'motif': [lead_tr_motif],
                },

            )
        dfs.append(results_df)
    final_df = pd.concat(dfs)
    final_df.to_csv('gs://cpg-bioheart-test/str/CD4_TCM_lead_TR_coordinates.csv', index = False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments


