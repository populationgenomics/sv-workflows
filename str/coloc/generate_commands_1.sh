#!/bin/bash

# Define the cell types
celltypes=(gdT B_intermediate ILC Plasmablast dnT ASDC cDC1 pDC NK_CD56bright MAIT B_memory CD4_CTL CD4_Proliferating CD8_Proliferating HSPC NK_Proliferating cDC2 CD16_Mono Treg CD14_Mono CD8_TCM CD4_TEM CD8_Naive CD4_TCM NK CD8_TEM CD4_Naive B_naive)
<<<<<<< Updated upstream
phenos=("haematocrit" "haemoglobin_concentration" "lymphocyte_count" "lymphocyte_percent" "mean_corpuscular_haemoglobin" "mean_corpuscular_haemoglobin_concentration" "mean_corpuscular_volume" "mean_platelet_volume" "mean_sphered_cell_volume" "neutrophil_count" "neutrophil_percent" "platelet_count" "platelet_crit" "platelet_distribution_width" "red_blood_cell_count" "red_blood_cell_distribution_width" "white_blood_cell_count")

# Iterate over the cell types and execute the command
for celltype in "${celltypes[@]}"; do
    for pheno in "${phenos[@]}"; do
        command="analysis-runner --dataset \"bioheart\" --description \"Run coloc for eGenes identified by STR analysis\" --access-level \"test\" --memory='8G' --image \"australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e\" --output-dir \"str/associatr\" coloc_ukbb_runner.py --pheno-output-name=gymrek-ukbb-$pheno --celltypes \"$celltype\" --max-parallel-jobs 10000"

        eval $command
    done
done
=======

# Define phenotypes with correct bash array syntax
phenos=("gymrek-ukbb-alanine-aminotransferase" "gymrek-ukbb-albumin" "gymrek-ukbb-alkaline_phosphatase" "gymrek-ukbb-apolipoprotein_a" "gymrek-ukbb-apolipoprotein_b" "gymrek-ukbb-aspartate_aminotransferase" "gymrek-ukbb-c_reactive_protein" "gymrek-ukbb-calcium" "gymrek-ukbb-cholesterol" "gymrek-ukbb-creatinine" "gymrek-ukbb-cystatin_c" "gymrek-ukbb-eosinophil_count" "gymrek-ukbb-eosinophil_percent" "gymrek-ukbb-gamma_glutamyltransferase" "gymrek-ukbb-glucose" "gymrek-ukbb-haematocrit" "gymrek-ukbb-haemoglobin_concentration" "gymrek-ukbb-lymphocyte_count" "gymrek-ukbb-lymphocyte_percent" "gymrek-ukbb-mean_corpuscular_haemoglobin" "gymrek-ukbb-mean_corpuscular_haemoglobin_concentration" "gymrek-ukbb-mean_corpuscular_volume" "gymrek-ukbb-mean_platelet_volume" "gymrek-ukbb-mean_sphered_cell_volume" "gymrek-ukbb-neutrophil_count" "gymrek-ukbb-neutrophil_percent" "gymrek-ukbb-platelet_count" "gymrek-ukbb-platelet_crit" "gymrek-ukbb-platelet_distribution_width" "gymrek-ukbb-red_blood_cell_count" "gymrek-ukbb-red_blood_cell_distribution_width" "gymrek-ukbb-white_blood_cell_count")

# Convert celltypes array to comma-separated string for the command
celltypes_str=$(IFS=,; echo "${celltypes[*]}")

# Iterate over the phenotypes and execute the command
for pheno in "${phenos[@]}"; do
        command="analysis-runner --dataset bioheart \
    --description Parse_coloc_results \
    --access-level test \
    --image australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e \
    --output-dir str/associatr \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-bioheart-test-analysis/str/associatr/coloc/sig_str_and_gwas_hit \
    --celltypes=$celltypes_str \
    --phenos=$pheno
"

        eval "$command"
done
>>>>>>> Stashed changes
