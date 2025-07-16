# TRTools Utilities for STR Workflows

This directory provides lightweight wrapper scripts around the [TRTools](https://github.com/gymrek-lab/TRTools) suite, enabling seamless integration of VCF preparation, merging, QC, summary statistics, and comparison into our Hail Batch + GCP STR pipelines.

---

## Contents

| Script                       | Purpose                                                                                              |
|------------------------------|------------------------------------------------------------------------------------------------------|
| **`merge_str_prep.py`**      | Reheaders and zips TR VCFs (e.g., ExpansionHunter output) for merging    |
| **`merge_str_runner.py`**    | Invoke TRTools’ `mergeSTR` to combine multiple prepped VCFs into one multi‐sample STR VCF            |
| **`qc_str_runner.py`**       | Run TRTools’ `qcSTR` to compute per‐sample and per‐locus QC metrics (call‐rate, stutter, etc.)        |
| **`stat_str_runner.py`**     | Run TRTools’ `statSTR` to generate per‐locus allele‐frequency and genotype‐distribution summaries     |
| **`compare_str_runner.py`**  | Run TRTools’ `compareSTR` to compare two STR callsets and report concordance metrics per‐locus/sample |


---



