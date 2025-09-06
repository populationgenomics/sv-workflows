# MCV Fine-Mapping Utilities

This directory implements the **Multiple-Causal-Variant (MCV)** fine-mapping framework in SuSIE using `susie()` for TR+SNV association signals. SuSIE estimates posterior inclusion probabilities (PIPs) and credible sets under a Bayesian model.

---

## Contents

| Script                           | Purpose                                                                                                    |
|----------------------------------|------------------------------------------------------------------------------------------------------------|
| `files_prep_dosages.py`                    | Obtains genotype dosages for the relevant cis windows                                  |
| `files_prep_residualizer.py`   | Regress out covariates from genotype dosages and pseudobulk expression values                  |
| `susie_mcv_runner.py`                  | Run `susie()` based on outputs from `files_prep_residualizer.py`                       |

---
