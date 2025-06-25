#!/usr/bin/env python3

"""
This script shards the UKBB SNP+STR catalog into chromosome-specific files.

analysis-runner --dataset "bioheart" --description "Liftover variants from hg19 to hg38" --access-level "test" \
    --output-dir "str/associatr/liftover" \
    --memory "4G" \
    --storage "5G" \
    catalog_sharder.py
"""

import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def sharder(phenotype):
    import gzip

    gwas_file = f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs_v4/white_british_{phenotype}_snp_str_gwas_results_hg38.tab.gz'
    with gzip.open(to_path(gwas_file), 'rb') as f:
        gwas = pd.read_csv(
            f,
            sep='\t',
        )

    for chrom_num in range(1, 23):
        chrom = f'chr{chrom_num}'
        gwas_chrom = gwas[gwas['chromosome'] == chrom]
        gwas_chrom.to_csv(
            f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs_v4/chr-specific/white_british_{phenotype}_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            sep='\t',
            index=False,
        )


def main():
    b = get_batch(name='catalog sharder')
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
        liftover_job = b.new_python_job('Catalog sharder for' + pheno)
        liftover_job.memory('32G')
        liftover_job.storage('10G')
        liftover_job.call(sharder, pheno)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
