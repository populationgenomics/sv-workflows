# associaTR Pipeline

This directory implements the full pseudobulk association pipeline for STR data, from aggregating single-cell counts through to multiple-testing correction. All steps are orchestrated via Hail Batch on GCP.

---

## Contents

| Script / Folder                      | Purpose                                                                                                  |
|--------------------------------------|----------------------------------------------------------------------------------------------------------|
| **`pseudobulk.py`**                  | Aggregate per-cell counts into per-individual, per-cell-type pseudobulk expression matrices.             |
| **`get_covariates.py`**                  | Extracts covariates (age, sex, genotype and RNA PCs).             |
| **`get_cis_numpy_files.py`**         | Extract per-gene phenotype & covariate NumPy arrays for cis-association testing.       |
| **`associatr_runner.py`**            | Main Hail-Batch driver to run associaTR per gene, per cell-type.                                         |
| **`meta_analysis/`**                 | Optional scripts to meta-analyze association results across cell types or cohorts.                       |
| **`multiple_testing_correction/`**   | Scripts to compute gene-level p-values (ACAT/Bonferroni) and FDR q-values (Storey).                     |
