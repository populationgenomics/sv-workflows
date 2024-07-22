#!/bin/bash

# Define the cell types
celltypes=(gdT B_intermediate ILC Plasmablast dnT ASDC cDC1 pDC NK_CD56bright MAIT B_memory CD4_CTL CD4_Proliferating CD8_Proliferating HSPC NK_Proliferating cDC2 CD16_Mono Treg CD14_Mono CD8_TCM CD4_TEM CD8_Naive CD4_TCM NK CD8_TEM CD4_Naive B_naive)

# Iterate over the cell types and execute the command
for celltype in "${celltypes[@]}"; do
    command="analysis-runner --dataset \"bioheart\" --description \"Run coloc for eGenes identified by STR analysis\" --access-level \"test\" --image \"australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e\" --output-dir \"str/associatr/fine_mapping/v3-ocv\" susie_runner.py  --celltypes \"$celltype\" --chromosomes "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22" --ld-dir "gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/prep_files/v2-whole-copies-only/correlation_matrix" --associatr-dir "gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results" --max-parallel-jobs 10000 --num-causal-variants=1 "

    eval $command
done
