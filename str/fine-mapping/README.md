# STR Fine-mapping Utilities

This directory contains scripts to perform fine-mapping of STR association signals using SuSiE, including pre-processing to clean indels and generate LD matrices.

---

## Contents

| Script                                                                                                                                     | Purpose                                                                                                  |
|--------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| [`remove_STR_indels.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/remove_STR_indels.py)               | Remove indels representing STRs and duplicate eSTR entries from associaTR outputs.                       |
| [`corr_matrix_maker.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/corr_matrix_maker.py)               | Compute variant correlation (LD) matrices for each locus from VCF or MatrixTable.                        |
| [`susie_runner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/susie_runner.py)                         | Run the SuSiE fine-mapping algorithm on summary statistics and LD matrices.                              |

---
