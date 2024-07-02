#!/usr/bin/env python3

"""
This script lifts over variants from hg19 to hg38.

The liftover is obtained from the Broad Institute liftover API.
It uses BCFTools liftover internally to get the liftover variant id.
New as of 2024 with support for multi-allelic variants.

analysis-runner --dataset "bioheart" --description "Liftover variants from hg19 to hg38" --access-level "test" \
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
    liftover_df = pd.read_csv(
        'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-gwas-catalogs/ukbb_snp_chr_pos_hg38_liftover.bed',
        sep='\t',
        header=None,
        names=['chromosome', 'position', 'end38', 'rsid'],
    )

    df = pd.merge(df, liftover_df, left_on='ID', right_on='rsid')
    df['varbeta'] = df['SE'] ** 2
    df['beta'] = df['BETA']
    df['snp'] = df['chromosome'] + '_' + df['position'].astype(str) + '_' + df['REF'] + '_' + df['ALT']
    df['p_value'] = df['P']

    # Write out the results
    df[['chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value']].to_csv(
        'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-gwas-catalogs/white_british_albumin_snp_gwas_results_hg38.tab.gz',
        sep='\t',
        index=False,
    )


def main():
    b = get_batch()
    liftover_job = b.new_python_job('Liftover variants from hg19 to hg38')
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
        liftover_job.memory('32G')
        liftover_job.storage('10G')
        liftover_job.call(liftover, pheno)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
