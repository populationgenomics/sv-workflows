# Collins et al. 2020 - A structural variation reference for med & pop genetics

## Data Summary

- Samples
  - 14,891 genomes (32x avg coverage)
    - 14,237 passed QC
      - Age median: 49 years
      - Ancestry: 46% Europeans, 35% Afroamericans
  - **12,653** after filtering out first-degree relatives
- SVs
  - Method: GATK-SV with four orthogonal evidence types
  - 433,371 SVs
    - Most low-quality SVs had incompletely resolved breakpoint junctions
  - **335,470** high-quality SVS
    - Median 7,439 SVs per genome
  - Most SVs were small (median 331 bp) and rare (92% of SVs had `AF < 1%`)
  - 50% of SVs were singletons (i.e. only one allele observed across all
    samples)
