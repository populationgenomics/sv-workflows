#!/bin/bash

# Define the cell types as a comma-separated string
#celltypes='B_intermediate,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,gdT,dnT,CD4_TCM'
celltypes='CD4_TEM,CD8_Naive'
# Convert the cell types into an array
IFS=',' read -ra celltype_array <<< "$celltypes"

# Loop through each cell type
for celltype in "${celltype_array[@]}"; do

    analysis-runner --dataset "bioheart" \
        --description "Run coloc for eGenes identified by STR analysis" \
        --access-level "test" \
        --memory "8G" \
        --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
        --output-dir "tenk10k/str/associatr/final_freeze" \
        coloc_runner.py \
        --snp-gwas-file=gs://cpg-bioheart-test/str/Trujillo_methylation_eQTLs/hg38_STRs_SNVs_parsed.tsv \
        --pheno-output-name="Trujillo_methylation_eQTLs" \
        --max-parallel-jobs=10000 \
        --celltypes="$celltype"
done
