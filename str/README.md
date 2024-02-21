# STR workflow

This folder contains scripts to genotype and characterise STRs using Hail Batch and GCP.
Runner scripts support HipSTR, GangSTR, and ExpansionHunter; however, downstream merging and QC scripts support ExpansionHunter VCFs only and are compatible with both sharded and unsharded catalogs.

## Suggested workflow for ExpansionHunter

- Genotype STRs using `runners/str_iterative_eh_runner.py`. If the catalog is >200k loci, we recommend sharding the catalog into 100k loci chunks.
- Prepare VCFs for multi-sample merging using `trtools/merge_str_prep.py`. If using a sharded catalog, this step should run on each sharded VCF of each sample in the cohort.
- Create a merged VCF using `trtools/merge_str_runner.py`. If using a sharded catalog, this step is performed on each shard separately.
- If using a sharded catalog, concatenate each sharded merged VCF using `helper/merge_str_vcf_combiner.py`.
- Final output is one VCF containing genotypes for all samples at all loci specified in the catalog, ready for import into Hail query:
  - BGZIP the VCF using `helper/bgzip.py`
  - Create a matrix table using `helper/mt_extractor.py`.

