#!/usr/bin/env bash

# Moves integrated anndata object from main-upload to main bucket

set -ex

gsutil -m cp gs://cpg-bioheart-main-upload/pseudobulk gs://cpg-bioheart-main/saige-qtl/integrated_anndata_object_from_HPC/224_libraries_concatenated_gene_info_donor_info.h5ad