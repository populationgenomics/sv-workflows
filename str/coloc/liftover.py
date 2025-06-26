#!/usr/bin/env python3

"""
This script lifts over variants from hg19 to hg38.

The liftover is obtained from the Broad Institute liftover API.
It uses BCFTools liftover internally to get the liftover variant id.
New as of 2024 with support for multi-allelic variants.

analysis-runner --dataset "tenk10k" --description "Liftover variants from hg19 to hg38" --access-level "test" \
    --output-dir "str/associatr/liftover" \
    --memory "4G" \
    --storage "5G" \
    liftover.py
"""

import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def liftover(phenotype):
    file_path = to_path(
        f'gs://cpg-bioheart-test-upload/str/ukbb-snp-catalogs/white_british_{phenotype}_snp_gwas_results.tab.gz',
    )
    df = pd.read_csv(file_path, sep='\t', compression='gzip', usecols=['ID', 'REF', 'ALT', 'BETA', 'SE', 'P'])
    df['composite_key'] = df['ID'] + '_' + df['REF'] + '_' + df['ALT'] # there are some multi-allelic variants with the same ID but different REF/ALT alleles

    columns = ['chr', 'position', 'refhg38', 'althg38', 'info']
    liftover_df = pd.read_csv(
        'gs://cpg-bioheart-test/str/gymrek_ukbb_snp_gwas_catalogs_v5/ukbb_snps_v2_lo_variants.tsv',
        sep='\t',
        header=None,
        names=columns
    )
    liftover_df['SRC_ID'] = liftover_df['info'].str.extract(r'SRC_ID=([^;]+)')
    liftover_df['SRC_REF_ALT'] = liftover_df['info'].str.extract(r'SRC_REF_ALT=([^;]+)')

    # Step 2: Split SRC_REF_ALT into SRC_REF and SRC_ALT, and assign them to new columns
    liftover_df[['SRC_REF', 'SRC_ALT']] = liftover_df['SRC_REF_ALT'].str.split(',', expand=True)
    liftover_df['composite_key'] = (
        liftover_df['SRC_ID'] + '_' + liftover_df['SRC_REF'] + '_' + liftover_df['SRC_ALT']
    )

    df = pd.merge(df, liftover_df, on='composite_key')  # noqa: PD015
    df['varbeta'] = df['SE'] ** 2
    df['beta'] = df['BETA']
    df['chromosome'] = df['chr']
    df['snp'] = df[['chr', 'position', 'refhg38', 'althg38']].astype(str).agg('_'.join, axis=1)
    df['p_value'] = df['P']

    #drop entries with duplicate SNPs (this occurs when 2 variants in hg19 map to the same coordinate in hg38 with the same REF/ALT alleles)
    df = df[~df['snp'].duplicated(keep=False)]


    # Write out the results
    df[['chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value']].to_csv(
        f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-gwas-catalogs_v6/white_british_{phenotype}_snp_gwas_results_hg38.tab.gz',
        sep='\t',
        index=False,
    )


def main():
    b = get_batch(name='liftover')
    phenotypes = [
        "alanine_aminotransferase",
        "albumin",
        "alkaline_phosphatase",
        "apolipoprotein_a",
        "apolipoprotein_b",
        "aspartate_aminotransferase",
        "c_reactive_protein",
        "calcium",
        "cholesterol",
        "creatinine",
        "cystatin_c",
        "eosinophil_count",
        "eosinophil_percent",
        "gamma_glutamyltransferase",
        "glucose",
        "glycated_haemoglobin",
        "haematocrit",
        "haemoglobin_concentration",
        "hdl_cholesterol",
        "igf_1",
        "ldl_cholesterol_direct",
        "lymphocyte_count",
        "lymphocyte_percent",
        "mean_corpuscular_haemoglobin",
        "mean_corpuscular_haemoglobin_concentration",
        "mean_corpuscular_volume",
        "mean_platelet_volume",
        "mean_sphered_cell_volume",
        "neutrophil_count",
        "neutrophil_percent",
        "phosphate",
        "platelet_count",
        "platelet_crit",
        "platelet_distribution_width",
        "red_blood_cell_count",
        "red_blood_cell_distribution_width",
        "shbg",
        "total_bilirubin",
        "total_protein",
        "triglycerides",
        "urate",
        "urea",
        "vitamin_d",
        "white_blood_cell_count",
    ]
    for pheno in phenotypes:
        liftover_job = b.new_python_job('Liftover variants from hg19 to hg38: ' + pheno)
        liftover_job.cpu(16)
        liftover_job.memory('highmem')
        liftover_job.storage('30G')
        liftover_job.call(liftover, pheno)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
