#!/bin/bash

# Define the cell types as a comma-separated string
celltypes='B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,gdT,dnT,CD4_TCM'
#celltypes='B_intermediate'

# Convert the cell types into an array
IFS=',' read -ra celltype_array <<< "$celltypes"

# Define the phenotypes as an array
pheno_names=(
    "alzheimer_GCST90027158"
    "breastca_GCST004988"
    "colorectalca_GCST90129505"
    "covid_GCST011071"
    "ibd_liu2023"
    "NHL_GCST90011819"
    "lungca_GCST004748"
    "lymphoma_GCST90018878"
    "parkinson_GCST009325"
    "prostateca_GCST90274713"
    "ra_GCST90132223"
    "sle_GCST003156"
    "myeloproliferative_GCST90000032"
    "lymphocytic_leukemia_GCST90011814"
    "nephrotic_GCST90258619"
    "kiryluk_IgAN"
)

# Loop through each cell type and phenotype
for pheno in "${pheno_names[@]}"; do
    analysis-runner --dataset "tenk10k" \
    --description "Parse coloc results" \
    --access-level "test" \
    --output-dir "str/associatr/final_freeze/meta_fixed/coloc-snp-only/sig_str_and_gwas_hit" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/meta_fixed/coloc-snp-only/sig_str_and_gwas_hit \
    --celltypes "$celltypes" \
    --phenos="$pheno"
done
