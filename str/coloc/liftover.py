#!/usr/bin/env python3

"""
This script lifts over variants from hg19 to hg38.

The liftover is obtained from the Broad Institute liftover API.
It uses BCFTools liftover internally to get the liftover variant id.
New as of 2024 with support for multi-allelic variants.

analysis-runner --dataset "bioheart" --description "Liftover variants from hg19 to hg38" --access-level "test" \
    --output-dir "str/associatr/liftover" \
    --memory "64G" \
    --storage "20G" \
    liftover.py --variants-file=gs://cpg-bioheart-test-upload/str/ukbb-snp-catalogs/white_british_cholesterol_snp_gwas_results.tab.gz
"""

import click
import pandas as pd
import logging
import requests
import time
from cpg_utils import to_path

@click.command()
@click.option('--variants-file', required=True)
def main(variants_file: str):
    """
    Opens the variants_file json and read the liftover variant id for each variant.
    If the liftover variant id is null, then get the liftover variant id
    from the Broad institiute API. Save the variants json with the new liftover ids.
    """
    file_path = to_path(variants_file)

    df = pd.read_csv(file_path, sep='\t', compression='gzip', usecols=['ID','REF','ALT','BETA','SE','P'])
    liftover_df = pd.read_csv('gs://cpg-bioheart-test/str/gymrek-ukbb-snp-gwas-catalogs/ukbb_snp_chr_pos_hg38_liftover.bed', sep = '\t', header=None, names=['chromosome', 'position', 'end38', 'rsid'])

    df = pd.merge(df, liftover_df, left_on = 'ID', right_on = 'rsid')
    df['varbeta'] = df['SE']**2
    df['beta'] = df['BETA']
    df['snp'] = df['chromosome'] + '_' + df['position'].astype(str) + '_' + df['REF'] + '_' + df['ALT']
    df['p_value'] = df['P']


    # Write out the results
    df[['chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value']].to_csv('gs://cpg-bioheart-test/str/gyremk-ukbb-snp-gwas-catalogs/white_british_cholesterol_snp_gwas_results_hg38.tab.gz', sep='\t', index=False)


if __name__ == '__main__':

    main()  # pylint: disable=no-value-for-parameter