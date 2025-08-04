#!/usr/bin/env python3
"""
This script performs pseudobulk (mean aggregation) of an input AnnData object (specifically by medium/low/high quality of a pathway)
Prior to pseudobulking, the following steps are performed:
- Filtering out lowly expressed genes
- Normalisation
- Log transformation (ln(1+x))
- Batch correction

Output is a TSV file by cell-type and chromosome-specific. Each row is a sample and each column is a gene.

analysis-runner --config pseudobulk.toml --access-level test --dataset tenk10k --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 --description "pseudobulk" --output-dir "str/cellstate/input_files/meanpool/stratified/tob/pseudobulk" python3 pseudobulk_macro.py

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


def pseudobulk(input_dir, chromosome,id_file_path, target_sum, min_pct,pathway):
    """
    Performs pseudobulking in a cell-type chromosome-specific manner
    """
    input_file_path_bint = f"{input_dir}/B_intermediate_chr{chromosome}.h5ad"
    adata_bint = sc.read_h5ad(to_path(input_file_path_bint))

    input_file_path_bmem = f"{input_dir}/B_memory_chr{chromosome}.h5ad"
    adata_bmem = sc.read_h5ad(to_path(input_file_path_bmem))

    adata = adata_bmem.concatenate(adata_bint, join="outer")

    # retain only samples in the id file
    with to_path(id_file_path).open() as csvfile:
        cpg_ids = list(csv.reader(csvfile))
    cpg_ids = cpg_ids[0]  # the function above creates a list in a list
    adata = adata[adata.obs['cpg_id'].isin(cpg_ids)]

    # read in adata pathway annotations
    pathway_annot = sc.read_h5ad(to_path('gs://cpg-bioheart-test/str/trdeepid/TR_b_cells_0408/meanpooling/pathway_attributions_macrosublabels.h5ad'))
    pathway_annot = pathway_annot[pathway_annot.obs["wg2_scpred_prediction"].isin(['B_intermediate', 'B_memory'])].copy() #filter to cell type
    subtype_series = pathway_annot.obs[pathway] #select the specific pathway's annotation
    adata.obs[pathway] = subtype_series #add the pathway annotation to the adata object

    # filter out lowly expressed genes
    n_all_cells = len(adata.obs.index)
    min_cells = math.ceil((n_all_cells * min_pct) / 100)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    for activity in ['low', 'medium', 'high']:
        # filter adata to only include cells with the current annotation
        adata_activity = adata[adata.obs[pathway] == activity].copy()

        # normalisation
        # we can't use pp.normalize_total because the expected input is a chr-specific anndata file, while
        # pp.normalize_total expects a whole genome-wide anndata file
        # obs.total_counts is the sum of counts (genome-wide) for each cell, determined in prior QC processing steps
        total_counts = adata_activity.obs['total_counts'].values
        scaled_x = (adata_activity.X / total_counts[:, np.newaxis]) * target_sum
        adata_activity.X = scaled_x

        # log transformation
        sc.pp.log1p(adata_activity)
        adata_activity.raw = adata_activity

        # batch correction
        sc.pp.regress_out(adata_activity, keys='sequencing_library')

        # initialise array to hold pseudobulk anndata objects
        pbs = []

        for sample in adata_activity.obs.individual.unique():
            samp_cell_subset = adata_activity[adata_activity.obs['individual'] == sample]
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

        chrom_num = input_file_path_bint.split('/')[-1].split('_chr')[1].split('.')[0]

        data_df.to_csv(
            output_path(f'B_int_mem/B_int_mem_chr{chrom_num}_{activity}_{pathway}_pseudobulk.csv'),
            index=False,
        )



def main():
    """
    Perform pseudobulk (mean aggregation) of an input AnnData object

    """
    b = get_batch('Run pseudobulk')

    for chromosome in get_config()['pseudobulk']['chromosomes'].split(','):
        j = b.new_python_job(name=f'Pseudobulk for B_int_mem: chr{chromosome}')
        j.image(image_path('scanpy'))
        j.cpu(get_config()['pseudobulk']['job_cpu'])
        j.memory(get_config()['pseudobulk']['job_memory'])
        j.storage(get_config()['pseudobulk']['job_storage'])
        j.call(
            pseudobulk,
            get_config()['pseudobulk']['input_dir'],
            chromosome,
            get_config()['pseudobulk']['sample_id_file_path'],
            get_config()['pseudobulk']['target_sum'],
            get_config()['pseudobulk']['min_pct'],
            get_config()['pseudobulk']['pathway']
        )
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
