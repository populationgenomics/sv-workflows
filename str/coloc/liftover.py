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
    columns = ['chr', 'pos', 'ref', 'alt', 'info']
    df = pd.read_csv(file_path, sep='\t', compression='gzip', usecols=['ID', 'REF', 'ALT', 'BETA', 'SE', 'P'])

    columns = ['chr', 'position', 'refhg38', 'althg38', 'info']
    liftover_df = pd.read_csv(
        'gs://cpg-bioheart-test/str/gymrek_ukbb_snp_gwas_catalogs_v2/ukbb_snps_hg38.sorted.chr_lo_variants.tsv',
        sep='\t',
        header=None,
        names=columns
    )
    liftover_df['SRC_ID'] = liftover_df['info'].str.extract(r'SRC_ID=([^;]+)')

    df = pd.merge(df, liftover_df, left_on='ID', right_on='SRC_ID')  # noqa: PD015
    df['varbeta'] = df['SE'] ** 2
    df['beta'] = df['BETA']
    df['chromosome'] = df['chr']
    df['snp'] = df['chr'] + '_' + df['position'].astype(str) + '_' + df['refhg38'] + '_' + df['althg38']
    df['p_value'] = df['P']

    # Write out the results
    df[['chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value']].to_csv(
        f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-gwas-catalogs_v2/white_british_{phenotype}_snp_gwas_results_hg38.tab.gz',
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
    for pheno in ["white_blood_cell_count"]:
        liftover_job = b.new_python_job('Liftover variants from hg19 to hg38: ' + pheno)
        liftover_job.memory('64G')
        liftover_job.storage('30G')
        liftover_job.call(liftover, pheno)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
