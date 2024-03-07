#!/usr/bin/env python3
"""
This script performs pseudobulk (mean aggregation) of an input AnnData object
Prior to pseudobulking, the following steps are performed:
- Normalisation
- Log transformation (ln(1+x))
- Batch correction

Output is a TSV file by cell-type and chromosome-specific. Each row is a sample and each column is a gene.

analysis-runner --access-level test --dataset bioheart --image australia-southeast1-docker.pkg.dev/cpg-common/images-dev/scanpy_sctransform:4.0 --description "pseudobulk" --output-dir "str/associatr/input_files" pseudobulk.py

"""

import scanpy as sc
from scipy.sparse import issparse

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

def pyScTransform(adata, ofile_path=None):
    """
    Function to call scTransform from Python
    """
    import rpy2.robjects as ro
    import anndata2ri

    ro.r('library(Seurat)')
    ro.r('library(scater)')
    anndata2ri.activate()

    sc.pp.filter_genes(adata, min_cells=5)

    if issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()

    for key in adata.layers:
        if issparse(adata.layers[key]):
            if not adata.layers[key].has_sorted_indices:
                adata.layers[key].sort_indices()

    ro.globalenv['adata'] = adata

    ro.r('seurat_obj = as.Seurat(adata, counts="X", data = NULL, assay = NULL)')
    print(ro.r('print(seurat_obj)'))

    ro.r('res <- SCTransform(object=seurat_obj, assay = "originalexp",vars.to.regress = c("pct_counts_mt","batch"),return.only.var.genes = FALSE, do.correct.umi = FALSE)')

    corrected_counts = ro.r('res@assays$SCT@counts').T

    adata.X = corrected_counts
    if output_file:
        adata.write_h5ad(str(ofile_path))


def main():
    """
    Perform pseudobulk (mean aggregation) of an input AnnData object

    """
    b = get_batch()
    sctransform_job = b.new_python_job('Run scTransform')
    expression_h5ad_path = to_path('gs://cpg-bioheart-test/str/anndata/saige-qtl/anndata_objects_from_HPC/ASDC_chr21.h5ad').copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path)

    sctransform_job.call(pyScTransform, adata, sctransform_job.ofile)
    sctransform_job.ofile.add_extension('.h5ad')
    b.write_output(sctransform_job.ofile, output_path('sctransform_chr21_ASDC.h5ad'))





if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
