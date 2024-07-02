#!/usr/bin/env python3
"""
This script lifts over variants from hg19 to hg38.

The liftover is obtained from the Broad Institute liftover API.
It uses BCFTools liftover internally to get the liftover variant id.
New as of 2024 with support for multi-allelic variants.

analysis-runner --dataset "bioheart" --description "Liftover variants from hg19 to hg38" --access-level "test" \
    --output-dir "str/associatr/liftover" \
    liftover.py --variants_file=gs://cpg-bioheart-test-upload/str/ukbb-snp-catalogs/white_british_cholesterol_snp_gwas_results.tab.gz
"""

import click
import pandas as pd
import logging
import requests
import time

from cpg_utils import to_path

from cloudpathlib import CloudPath

LIFTOVER_URL_BASE = 'https://liftover-xwkwwwxdwq-uc.a.run.app/liftover/?hg={hg}&format=variant&chrom={chrom}&pos={pos}&ref={ref}&alt={alt}'


def get_broad_liftover(variant, hg):
    """
    Get Broad Institute liftover data for a variant.

    variant (required) a variant in the format "chrom-pos-ref-alt"
    hg (required) can be 37 or 38
    """
    variant = variant.split('-')
    assert len(variant) == 4

    if hg in ('19', '37', 'GRCh37'):
        hg = 'hg19-to-hg38'
    elif hg in ('38', 'GRCh38'):
        hg = 'hg38-to-hg19'
    else:
        raise ValueError('Reference genome must be 19/37 or 38')

    chrom, pos, ref, alt = variant

    url = LIFTOVER_URL_BASE.format(hg=hg, chrom=chrom, pos=pos, ref=ref, alt=alt)
    logging.info('Getting liftover data from %s', url)

    response = requests.get(url, verify=True)
    response.raise_for_status()

    data = response.json()
    liftover_pos = data['output_pos']
    liftover_ref = data['output_ref']
    liftover_alt = data['output_alt']

    return f'{chrom}-{liftover_pos}-{liftover_ref}-{liftover_alt}'


@click.command()
@click.option('--variants-file', required=True)
def main(variants_file: str):
    """
    Opens the variants_file json and read the liftover variant id for each variant.
    If the liftover variant id is null, then get the liftover variant id
    from the Broad institiute API. Save the variants json with the new liftover ids.
    """
    file_path = to_path(variants_file)

    df = pd.read_csv(file_path, sep='\t', compression='gzip')

    # Get the Broad liftover variant ID for each variant
    for row in df.iterrows():
        liftover_input = row['#CHROM'] + '-' + row['POS'] + '-' + row['REF'] + '-' + row['ALT']

        broad_liftover_hg38 = get_broad_liftover(liftover_input, '19')
        time.sleep(6)  # 12 requests per minute max
        row['position'] = broad_liftover_hg38.split('-')[1]
        row['chromosome'] = 'chr'+row['#CHROM'].astype(str)
        row['varbeta'] = row['SE']**2
        row['beta'] = row['BETA']
        row['snp'] = row['#CHROM'] + '_' + row['POS'].astype(str) + '_' + row['REF'] + '_' + row['ALT']
        row['p_value'] = row['P']

    row[['chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value']].to_csv('gs://cpg-bioheart-test/str/gyremk-ukbb-snp-gwas-catalogs/white_british_cholesterol_snp_gwas_results_hg38.tab.gz', sep='\t', index=False)


if __name__ == '__main__':

    main()  # pylint: disable=no-value-for-parameter