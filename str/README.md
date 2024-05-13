# STR workflow

This folder contains scripts to genotype and characterise STRs using Hail Batch and GCP.
Runner scripts support HipSTR, GangSTR, and ExpansionHunter; however, downstream merging and QC scripts support ExpansionHunter VCFs only and are compatible with both sharded and unsharded catalogs.

## ExpansionHunter genome-wide genotyping

- Genotype STRs using `runners/str_iterative_eh_runner.py`. If the catalog is >200k loci, we recommend sharding the catalog into 100k loci chunks.
- Prepare VCFs for multi-sample merging using `trtools/merge_str_prep.py`. If using a sharded catalog, this step should run on each sharded VCF of each sample in the cohort.
- Create a merged VCF using `trtools/merge_str_runner.py`. If using a sharded catalog, this step is performed on each shard separately.
- If using a sharded catalog, concatenate each sharded merged VCF using `helper/merge_str_vcf_combiner.py`.
- Final output is one VCF containing genotypes for all samples at all loci specified in the catalog, ready for import into Hail query:
  - BGZIP the VCF using `helper/bgzip.py`
  - Create a matrix table using `helper/mt_extractor.py`.

## Pseudobulk scRNA association pipeline using associaTR

Assumes scRNA raw data have been processed and cells have been typed using the Powell lab pipeline, producing chromosome and cell-type specific h5ad objects.

- Perform pseudobulking (mean aggregation) of scRNA data using `pseudobulk.py`.
- Prepare necessary inputs for associaTR: 1) covariates using `get_covariates.py` and 2) numpy objects containing the covariates and phenotypes (pseudobulked expression) using `get_cis_numpy_files.py`.
- Annotate and apply QC filters to the STR genotypes matrix table using `qc/qc_annotator.py` and subsequently `qc/qc_filters_associatr.py`. The latter produces chromosome-specific VCFs for input into associaTR.
- Run associaTR with `associatr_runner.py`
- Optional: Meta-analysis runner scripts in `associatr/meta_analysis`
- Perform multiple testing correction (`associatr/multiple_testing_correction`)
  - at the gene-level using ACAT correction with `run_gene_level_pval.py`. Option to use Bonferroni correction instead.
  - control for FDR using Storey q-values with `run_storey.py`.
 
## Downstream analysis

### Colocalisation

- Runner scripts stored in `coloc`: 1) `coloc_runner_{phenotype}.py` and 2) `coloc_results_parser.py` which consolidates the outputs from 1) into one file.
- If colocalisation was done using SNP data, then we check whether the lead eSTR associated with the colocalized locus is in LD with at least one SNP in the GWAS catalog using 1) `ld_runner.py` and 2) `ld_results_parser.py` which consolidates the outputs from 1) into one file.
