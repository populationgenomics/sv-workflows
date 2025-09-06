# STR workflow

This folder contains scripts to genotype and characterise STRs using Hail Batch and GCP.  
Runner scripts support HipSTR, GangSTR, and ExpansionHunter; however, downstream merging and QC scripts support ExpansionHunter VCFs only and are compatible with both sharded and unsharded catalogs.  
Sample VCFs can be obtained from the [TRTools Git](https://github.com/gymrek-lab/TRTools/tree/master/example-files) and are fully compatible with this pipeline. Users may run the workflow with these files to verify installation and become familiar with the expected input format before applying it to their own data.

---

## ExpansionHunter genome-wide genotyping

- Genotype STRs using [`runners/str_iterative_eh_runner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/runners/str_iterative_eh_runner.py).  
  If the catalog is >200k loci, we recommend sharding the catalog into 100k loci chunks.
- Prepare VCFs for multi-sample merging using [`trtools/merge_str_prep.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/trtools/merge_str_prep.py).  
  If using a sharded catalog, this step should run on each sharded VCF of each sample in the cohort.
- Create a merged VCF using [`trtools/merge_str_runner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/trtools/merge_str_runner.py).  
  If using a sharded catalog, this step is performed on each shard separately.
- If using a sharded catalog, concatenate each sharded merged VCF using [`helper/merge_str_vcf_combiner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/helper/merge_str_vcf_combiner.py).
- Final output is one VCF containing genotypes for all samples at all loci specified in the catalog, ready for import into Hail:
  - BGZIP the VCF using [`helper/bgzip.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/helper/bgzip.py)
  - Create a matrix table using [`helper/mt_extractor.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/helper/mt_extractor.py).

---

## Pseudobulk scRNA association pipeline using associaTR

Assumes scRNA raw data have been processed and cells have been typed using the [Powell lab pipeline](https://github.com/powellgenomicslab/tenk10k_phase1), producing chromosome and cell-type specific h5ad objects.

- Perform pseudobulking (mean aggregation) of scRNA data using [`pseudobulk.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/pseudobulk.py).
- Prepare necessary inputs for associaTR:
  - covariates using [`get_covariates.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/get_covariates.py)
  - numpy objects containing the covariates and phenotypes (pseudobulked expression) using [`get_cis_numpy_files.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/get_cis_numpy_files.py).
- **Note:** Conditional analysis (conditioning on the lead STR or SNP signal) is also specified in [`get_cis_numpy_files.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/get_cis_numpy_files.py).
- Annotate and apply QC filters to the STR genotypes matrix table using:
  - [`qc/qc_annotator.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/qc/qc_annotator.py)
  - [`qc/qc_filters_associatr.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/qc/qc_filters_associatr.py).  
    This produces chromosome-specific VCFs for input into associaTR.
- Run associaTR with [`associatr_runner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/associatr_runner.py).
- Optional: Meta-analysis runner scripts in [`associatr/meta_analysis`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/meta_analysis).
- Perform multiple testing correction (`associatr/multiple_testing_correction`):
  - Gene-level ACAT correction with [`run_gene_level_pval.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/multiple_testing_correction/run_gene_level_pval.py) (option to use Bonferroni).
  - FDR control using Storey q-values with [`run_storey.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/multiple_testing_correction/run_storey.py).

---

## Running SNPs in associaTR

To run SNPs using the associaTR pipeline:

1. Subset the SNP VCF to target samples using [`vcf_sample_subsetter.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/helper/vcf_sample_subsetter.py)
2. Format the SNP VCF into an EH-style VCF using [`snp_vcf_for_associatr.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/helper/snp_vcf_for_associatr.py)
3. Bgzip and tabix the VCFs with [`bgzip_tabix.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/helper/bgzip_tabix.py)
4. Use these input VCF files with [`associatr_runner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/associatr_runner.py).

Finally, combine eSTR and eSNP associaTR results using [`dataframe_concatenator.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/associatr/meta_analysis/dataframe_concatenator.py).

---

## Downstream analysis

### MashR
- To do once scripts in main 

### Cell specificity 
- To do once scripts in main 

### Finemapping (`fine-mapping`) with SuSIE (multiple causal variant assumption)

- Run [`remove_STR_indels.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/remove_STR_indels.py) to remove indels that represent STRs, as well as duplicate eSTRs in the associaTR outputs.
- Use [`files_prep_dosages.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/mcv/files_prep_dosages.py)) to obtain the genotypes from the VCFs.
- Run [`files_prep_residualizer.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/mcv/files_prep_residualizer.py) to regress out the covariates from the genotypes and the pseudobulked exprssion values. 
- Run SusieR with [`susie_mcv_runner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/fine-mapping/mcv/susie_mcv_runner.py).

### Colocalisation (`coloc`)

- Various runner scripts stored in [`coloc/`](https://github.com/populationgenomics/sv-workflows/blob/main/str/coloc/) (format: `coloc_runner_{phenotype}.py`).
- Consolidate outputs with [`coloc_results_parser.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/coloc/coloc_results_parser.py).
- If colocalisation was performed using SNP data, check LD with lead eSTRs using [`coloc_ld_runner.py`](https://github.com/populationgenomics/sv-workflows/blob/main/str/coloc/coloc_ld_runner.py).

---

📌 **Note:** If you encounter issues or have questions, please open an issue on this repository.  

## Citation and data availability

- Preprint discussing the work in detail can be accessed at doi.org/10.1101/2024.11.02.621562. 

 - Raw summary association statistics and the polymorphic variant catalog will be available at 10.5281/zenodo.15009519 upon formal publication. 

