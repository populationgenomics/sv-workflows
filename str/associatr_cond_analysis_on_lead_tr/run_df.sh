#!/bin/bash

# Define an array of cell types
celltypes=(
  "B_intermediate" "Plasmablast" "ASDC" "cDC1" "pDC" "NK_CD56bright" "MAIT" "B_memory" "CD4_CTL"
  "CD4_Proliferating" "CD8_Proliferating" "HSPC" "NK_Proliferating" "cDC2" "CD16_Mono" "Treg"
  "CD14_Mono" "CD8_TCM" "CD4_TEM" "CD8_Naive" "NK" "CD8_TEM" "CD4_Naive" "B_naive" "gdT" "dnT"
  "CD4_TCM"
)
# Iterate over each cell type and run the analysis
for celltype in "${celltypes[@]}"; do
  analysis-runner --dataset "bioheart" --description "concatenate meta-analysis results" --access-level "test" \
    --output-dir "tenk10k/str/associatr/final_freeze/cond_analysis_on_tr/snps_and_strs/bioheart_n975_and_tob_n950/meta_results" \
    dataframe_concatenator.py \
    --input-dir-1=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/cond_analysis_on_tr/common_variants_snps/bioheart_n975_and_tob_n950/meta_results \
    --input-dir-2=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/cond_analysis_on_tr/bioheart_n975_and_tob_n950/meta_results \
    --celltypes="$celltype"
done