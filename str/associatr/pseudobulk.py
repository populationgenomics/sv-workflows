#!/usr/bin/env python3
"""
This script performs pseudobulk (mean aggregation) of an input AnnData object
Prior to pseudobulking, the following steps are performed:
- Filtering out lowly expressed genes
- Normalisation
- Log transformation (ln(1+x))
- Batch correction

Output is a TSV file by cell-type and chromosome-specific. Each row is a sample and each column is a gene.

analysis-runner --config pseudobulk.toml --access-level test --dataset bioheart --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 --description "pseudobulk" --output-dir "str/test" python3 pseudobulk.py

"""
import csv
import logging
import math

import numpy as np
import pandas as pd
import scanpy as sc

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, image_path, output_path


def pseudobulk(input_file_path, id_file_path, target_sum, min_pct):
    """
    Performs pseudobulking in a cell-type chromosome-specific manner
    """
    expression_h5ad_path = to_path(input_file_path).copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path)

    # retain only samples in the id file
    with to_path(id_file_path).open() as csvfile:
        cpg_ids = list(csv.reader(csvfile))
    cpg_ids = cpg_ids[0]  # the function above creates a list in a list
    adata = adata[adata.obs['cpg_id'].isin(cpg_ids)]

    # filter out lowly expressed genes
    n_all_cells = len(adata.obs.index)
    min_cells = math.ceil((n_all_cells * min_pct) / 100)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # normalisation
    # we can't use pp.normalize_total because the expected input is a chr-specific anndata file, while
    # pp.normalize_total expects a whole genome-wide anndata file
    # obs.total_counts is the sum of counts (genome-wide) for each cell, determined in prior QC processing steps
    total_counts = adata.obs['total_counts'].values
    scaled_x = (adata.X / total_counts[:, np.newaxis]) * target_sum
    adata.X = scaled_x

    # log transformation
    sc.pp.log1p(adata)
    adata.raw = adata

    # batch correction
    sc.pp.regress_out(adata, keys='sequencing_library')

    # initialise array to hold pseudobulk anndata objects
    pbs = []

    for sample in adata.obs.individual.unique():
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

    cell_type = input_file_path.split('/')[-1].split('_chr')[0]
    chrom_num = input_file_path.split('/')[-1].split('_chr')[1].split('.')[0]

    data_df.to_csv(
        output_path(f'{cell_type}/{cell_type}_chr{chrom_num}_pseudobulk.csv'),
        index=False,
    )


# inputs:


def main():
    """
    Perform pseudobulk (mean aggregation) of an input AnnData object

    """
    b = get_batch('Run pseudobulk')

    for cell_type in get_config()['pseudobulk']['cell_types'].split(','):
        for chromosome in get_config()['pseudobulk']['chromosomes'].split(','):
            input_file = f"{get_config()['pseudobulk']['input_dir']}/{cell_type}_chr{chromosome}.h5ad"
            j = b.new_python_job(name=f'Pseudobulk for {cell_type}: chr{chromosome}')
            j.image(image_path('scanpy'))
            j.cpu(get_config()['pseudobulk']['job_cpu'])
            j.memory(get_config()['pseudobulk']['job_memory'])
            j.storage(get_config()['pseudobulk']['job_storage'])
            j.call(
                pseudobulk,
                input_file,
                get_config()['pseudobulk']['sample_id_file_path'],
                get_config()['pseudobulk']['target_sum'],
                get_config()['pseudobulk']['min_pct'],
            )
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
