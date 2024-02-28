#!/usr/bin/env python3
"""
This script performs pseudobulk (mean aggregation) of an input AnnData object
Prior to pseudobulking, the following steps are performed:
- Remove cells with no counts
- Normalisation
- Log transformation (ln(1+x))
- Batch correction using ComBat

Output is a TSV file by cell-type. Each row is a sample and each column is a gene.

analysis-runner --access-level test --dataset bioheart --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 --description "pseudobulk" --storage=200G --cpu=16 --memory=highmem --output-dir "str/associatr/input_files" pseudobulk.py --input-file=gs://cpg-bioheart-test/str/anndata/test.h5ad

"""
import click
import numpy as np
import pandas as pd
import scanpy as sc

from cpg_utils.hail_batch import output_path
from cpg_utils import to_path


# inputs:
@click.option('--input-file', help='GCS Path to the input AnnData object')
@click.command()
def main(input_file):
    """
    Perform pseudobulk (mean aggregation) of an input AnnData object

    Args
        input_file (str): GCS Path to the input AnnData object
    """
    expression_h5ad_path = to_path(input_file).copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path)

    # remove cells with no counts
    sc.pp.filter_cells(adata, min_counts=1)

    # normalisation
    sc.pp.normalize_total(adata, target_sum=1e6)
    # log transformation
    sc.pp.log1p(adata)
    adata.raw = adata

    # batch correction
    sc.pp.combat(adata, key='sequencing_library')

    for cell_type in adata.obs.wg2_scpred_prediction.unique():
        cell_type_adata = adata[adata.obs['wg2_scpred_prediction'] == cell_type]
        # each cell type will have its own pseudobulk object
        pbs = []
        for chr in cell_type_adata.var.chr.unique():
            chr_adata = cell_type_adata[:, cell_type_adata.var['chr'] == chr]
            for sample in chr_adata.obs.individual.unique():
                samp_cell_subset = adata[adata.obs['individual'] == sample]
                # pseudobulk mean aggregation
                rep_adata = sc.AnnData(
                    X=np.matrix(samp_cell_subset.X.mean(axis=0)),
                    var=samp_cell_subset.var[[]],
                )
                rep_adata.obs_names = [sample]
                pbs.append(rep_adata)
        # concatenate individual-specific pseudobulk objects
        pb = sc.concat(pbs)

        # convert pb into dataframe
        data_df = pd.DataFrame(pb.X, index=pb.obs.index, columns=pb.var.index)
        data_df = data_df.rename_axis('individual').reset_index()
        cell_type_name = cell_type.replace(' ', '_')
        data_df.to_csv(
            output_path(f'pseudobulk/{cell_type_name}_pseudobulk.csv'), index=False
        )


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
