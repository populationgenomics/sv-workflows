# STR QC Utilities

This directory contains scripts to **annotate** and **filter** STR genotype data for downstream association testing and analysis in our Hail Batch + GCP pipeline. These wrappers leverage the [TRTools](https://github.com/gymrek-lab/TRTools) QC functionality and custom filtering logic tailored for the associaTR workflow.

---

## Contents

| Script                          | Purpose                                                                                             |
|---------------------------------|-----------------------------------------------------------------------------------------------------|
| **`qc_annotator.py`**           | Annotate a Hail MatrixTable of STR genotypes with per-variant and per-sample QC metrics              |
| **`qc_filters_associatr.py`**   | Apply QC filters to the annotated MatrixTable and export chromosome-specific VCFs for associaTR use |

---

## `qc_annotator.py`

**Description**  
Annotates an ExpansionHunter MatrixTable (`.mt`) with allele- and locus-level repeat statistics.

**Annotations added**  
- `rep_length_alleles` Repeat length of each allele  
- `motif_length` Length of the STR motif  
- `bp_length_alleles` Base-pair length of each allele  
- `allele_1_rep_length`, `allele_2_rep_length` Repeat lengths of allele 1 & 2  
- `allele_1_bp_length`, `allele_2_bp_length` Base-pair lengths of allele 1 & 2  
- `aggregated_info.mode_allele` Mode allele at each locus  
- `num_alleles` Number of distinct alleles observed per locus  
- `allele_1_minus_mode`, `allele_2_minus_mode` Difference from mode allele  
- `sum_alleles_is_not_mode` Count of non-mode alleles at each locus  
- `prop_alleles_is_not_mode` Proportion of alleles that differ from the mode  
- `binom_hwep` Binomial HWE P-value per locus  
- `obs_het` Observed heterozygosity (proportion of het calls)

These annotations are stored in the MatrixTable’s row and column fields, ready for filtering.

## `qc_filters_associatr.py`

Filters a QC‐annotated ExpansionHunter MatrixTable and exports per‐chromosome VCFs ready for input into the `associatr_runner.py` pipeline.

## Description

This script takes the output of `qc_annotator.py` (a Hail MatrixTable with QC annotations) and applies stringent locus‐ and sample‐level filters before exporting chromosome‐specific VCF files. All filtering is performed in Hail for speed and scalability.

### Filters Applied

1. **Remove monomorphic variants**  
2. **Locus call rate ≥ 0.90**  
3. **Sample call rate ≥ 0.99**  
4. **Observed heterozygosity ≥ 0.00995** (≈ MAF 0.5%)  
5. **Hardy–Weinberg equilibrium P ≥ 1×10⁻⁶** (binomial HWE test)  
6. **Genotype range mask**: any allele call outside `[mode – 30, mode + 20]` is set to missing  
7. **Re‐enforce locus call rate ≥ 0.90** after masking
