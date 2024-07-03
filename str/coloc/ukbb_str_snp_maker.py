#!/usr/bin/env python3

"""
This script combines the STR an SNP Gymrek UKBB catalogs.
STR UKBB catalogs are parsed with a harmonised mapping file to ensure its definition matches the eQTL callset.

analysis-runner --dataset "bioheart" --description "Liftover variants from hg19 to hg38" --access-level "test" \
    --output-dir "str/associatr/liftover" \
    --memory "4G" \
    --storage "5G" \
    ukbb_str_snp_maker.py
"""

import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def liftover(phenotype):
    import gzip

    str_gwas_file = 'gs://cpg-bioheart-test/str/gymrek-ukbb-str-gwas-catalogs/gymrek-ukbb-str-gwas-catalogs/white_british_albumin_str_gwas_results.tab.gz'
    with gzip.open(to_path(str_gwas_file), 'rb') as f:
        str_gwas = pd.read_csv(
            f,
            sep='\t',
            usecols=[
                'chromosome',
                'beta',
                'standard_error',
                'p_value',
                'repeat_unit',
                'start_pos (hg38)',
                'end_pos (hg38)',
            ],
        )
    str_gwas['chromosome'] = 'chr' + str_gwas['chromosome'].astype(str)

    mapping_df = pd.read_csv('gs://cpg-bioheart-test/str/gymrek-ukbb-str-gwas-catalogs/ukbb_str_harmonised_mapping.csv')

    df = pd.merge(
        str_gwas,
        mapping_df,
        left_on=['chromosome', 'start_pos (hg38)', 'repeat_unit'],
        right_on=['chrom', 'gwas_pos', 'gwas_motif'],
    )

    df['varbeta'] = df['standard_error'] ** 2
    df['position'] = df['catalog_pos']

    # Write out the results
    df = df[['chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value']]

    snp_gwas_parsed_file = f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-gwas-catalogs/white_british_{phenotype}_snp_gwas_results_hg38.tab.gz'
    with gzip.open(to_path(snp_gwas_parsed_file), 'rb') as f_snp:
        snp_gwas = pd.read_csv(f_snp, sep='\t')

    df_combined = pd.concat([df, snp_gwas])
    df_combined.to_csv(
        f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/white_british_{phenotype}_snp_str_gwas_results_hg38.tab.gz',
        sep='\t',
        index=False,
    )


def main():
    b = get_batch(name='liftover')
    phenotypes = 'albumin'
    for pheno in phenotypes:
        liftover_job = b.new_python_job('Parse STR UKBB and combine with SNP ' + pheno)
        liftover_job.memory('64G')
        liftover_job.storage('10G')
        liftover_job.call(liftover, pheno)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
