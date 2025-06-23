#!/bin/bash

# Define the cell types as a comma-separated string
celltypes='B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,gdT,dnT,CD4_TCM'
#celltypes='B_intermediate'

# Convert the cell types into an array
IFS=',' read -ra celltype_array <<< "$celltypes"

# Define the phenotypes as an array
pheno_names=(
    "gymrek-ukbb-alanine_aminotransferase"
    "gymrek-ukbb-albumin"
    "gymrek-ukbb-alkaline_phosphatase"
    "gymrek-ukbb-apolipoprotein_a"
    "gymrek-ukbb-apolipoprotein_b"
    "gymrek-ukbb-aspartate_aminotransferase"
    "gymrek-ukbb-c_reactive_protein"
    "gymrek-ukbb-calcium"
    "gymrek-ukbb-cholesterol"
    "gymrek-ukbb-creatinine"
    "gymrek-ukbb-cystatin_c"
    "gymrek-ukbb-eosinophil_count"
    "gymrek-ukbb-eosinophil_percent"
    "gymrek-ukbb-gamma_glutamyltransferase"
    "gymrek-ukbb-glucose"
    "gymrek-ukbb-glycated_haemoglobin"
    "gymrek-ukbb-haematocrit"
    "gymrek-ukbb-haemoglobin_concentration"
    "gymrek-ukbb-hdl_cholesterol"
    "gymrek-ukbb-igf_1"
    "gymrek-ukbb-ldl_cholesterol_direct"
    "gymrek-ukbb-phosphate"
    "gymrek-ukbb-shbg"
    "gymrek-ukbb-total_bilirubin"
    "gymrek-ukbb-total_protein"
    "gymrek-ukbb-triglycerides"
    "gymrek-ukbb-urate"
    "gymrek-ukbb-urea"
    "gymrek-ukbb-vitamin_d"
    "gymrek-ukbb-white_blood_cell_count"
    "gymrek-ukbb-lymphocyte_count"
    "gymrek-ukbb-lymphocyte_percent"
    "gymrek-ukbb-mean_corpuscular_haemoglobin"
    "gymrek-ukbb-mean_corpuscular_haemoglobin_concentration"
    "gymrek-ukbb-mean_corpuscular_volume"
    "gymrek-ukbb-mean_platelet_volume"
    "gymrek-ukbb-mean_sphered_cell_volume"
    "gymrek-ukbb-neutrophil_count"
    "gymrek-ukbb-neutrophil_percent"
    "gymrek-ukbb-platelet_count"
    "gymrek-ukbb-platelet_crit"
    "gymrek-ukbb-platelet_distribution_width"
    "gymrek-ukbb-red_blood_cell_count"
    "gymrek-ukbb-red_blood_cell_distribution_width"
)

# Loop through each cell type and phenotype
for pheno in "${pheno_names[@]}"; do
    analysis-runner --dataset "tenk10k" \
    --description "Parse coloc results" \
    --access-level "test" \
    --output-dir "str/associatr/final_freeze/meta_fixed/coloc/sig_str_and_gwas_hit" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/meta_fixed/coloc/sig_str_and_gwas_hit \
    --celltypes "$celltypes" \
    --phenos="$pheno"
done
