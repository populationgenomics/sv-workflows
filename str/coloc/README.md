# STR Colocalisation Utilities

This directory contains scripts to perform colocalisation analyses between STR association signals (e.g. eSTRs) and external phenotype summary statistics (GWAS) using the **coloc** framework via Hail Batch.

---

## Contents

| Script                                      | Description                                                                                                                      |
|---------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------|
| `coloc_runner.py`               | Runner scripts to perform colocalisation for a given phenotype. The script launches per-chromosome or per-locus coloc analyses based on an input GWAS catalog. |
| `coloc_ukbb_runner.py`               | Runner scripts to perform colocalisation for a given UKBB phenotype. The script launches per-chromosome or per-locus coloc analyses based on an input GWAS catalog. |
| `coloc_results_parser.py`                   | Gathers all of the per-phenotype/run output files and merges them into a single summary table of colocalisation statistics (PP.H0…PP.H4, etc.).             |
| `coloc_ld_runner.py`                        | For colocalisations based on SNP data, computes LD between the lead STR and GWAS catalog variants at each locus to support interpretation.                      |
