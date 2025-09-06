# Runner Scripts for STR Genotyping

This directory contains **runner scripts for executing GangSTR, HipSTR, ExpansionHunter, and REViewer** on large-scale cohorts using [Hail Batch](https://hail.is/docs/batch/) and Google Cloud Platform (GCP).

These scripts are designed for scalability and reproducibility in distributed environments, enabling efficient processing of thousands of whole-genome samples.

---

## Contents

| Script                                   | Description                                                       |
|------------------------------------------|-------------------------------------------------------------------|
| `str_iterative_eh_runner.py`             | Run [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) iteratively across a catalog of STR loci, supporting catalog sharding for large-scale genotyping. |
| `str_iterative_gangstr_runner.py`                  | Run [GangSTR](https://github.com/gymreklab/GangSTR) on a set of samples with pre-defined STR regions. |
| `str_iterative_hipstr_runner.py`                   | Run [HipSTR](https://github.com/tfwillems/HipSTR) for multi-sample joint STR genotyping. |
| `reviewer_runner.py`                 | Run [REViewer](https://github.com/Illumina/REViewer) to visualize read-level support for repeat expansions from ExpansionHunter calls. |

---

## Features

- Built on **Hail Batch**, leveraging the scalability and fault tolerance of GCP.
- Containerized execution (Docker) to ensure reproducibility across environments.
- Support for sharded catalogs (e.g., for large ExpansionHunter runs on >200k loci).
- Standardized output paths compatible with downstream merging, QC, and association pipelines.

---
