#!/bin/bash

# Define the cell types
celltypes=(gdT B_intermediate ILC Plasmablast dnT ASDC cDC1 pDC NK_CD56bright MAIT B_memory CD4_CTL CD4_Proliferating CD8_Proliferating HSPC NK_Proliferating cDC2 CD16_Mono Treg CD14_Mono CD8_TCM CD4_TEM CD8_Naive CD4_TCM NK CD8_TEM CD4_Naive B_naive)

# Iterate over the cell types and execute the command
for celltype in "${celltypes[@]}"; do
    command="analysis-runner --dataset \"bioheart\" --description \"Run coloc for eGenes identified by STR analysis\" --access-level \"test\" --memory='8G' --image \"australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e\" --output-dir \"str/associatr\" coloc_ukbb_runner.py --pheno-output-name=gymrek-ukbb-alkaline_phosphatase --celltypes \"$celltype\" --max-parallel-jobs 10000"

    eval $command
done