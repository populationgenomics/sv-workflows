##cell_types=("ILC" "Plasmablast" "ASDC" "cDC1" "pDC" "NK_CD56bright" "MAIT" "B_memory" "CD4_CTL" "CD4_Proliferating" "CD8_Proliferating" "HSPC" "NK_Proliferating" "cDC2" "CD16_Mono" "Treg" "CD14_Mono" "CD8_TCM" "CD4_TEM" "CD8_Naive" "NK" "CD8_TEM" "CD4_Naive" "B_naive" "gdT" "dnT" "CD4_TCM" "B_intermediate" "CD4_TCM_permuted")
#cell_types=("CD4_TCM_permuted")
#"CD4_TCM_permuted" "B_intermediate"
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22")
#"chr6"
for chromosome in "${chromosomes[@]}"; do
analysis-runner --dataset "bioheart" \
    --description "Calculate LD between STR and SNPs" \
    --access-level "test" \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "tenk10k/str/associatr/final_freeze/fine_mapping/prep_files/v1" \
    corr_matrix_maker.py --snp-vcf-dir=gs://cpg-bioheart-test/tenk10k/str/associatr/common_variant_snps \
    --str-vcf-dir=gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific \
    --celltypes=gdT,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive,B_intermediate \
    --job-storage=30G --job-cpu=2 \
    --max-parallel-jobs=500 \
    --str-fdr-dir=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/meta_results/fdr_qvals/using_acat \
    --associatr-dir=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/snps_and_strs/bioheart_n975_and_tob_n950/rm_str_indels_dup_strs/meta_results \
    --chromosomes="$chromosome"
done