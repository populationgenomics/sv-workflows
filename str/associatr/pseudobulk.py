#!/usr/bin/env python3
"""
This script performs pseudobulk (mean aggregation) of an input AnnData object
Prior to pseudobulking, the following steps are performed:
- Filtering out lowly expressed genes
- Normalisation
- Log transformation (ln(1+x))
- Batch correction

Output is a TSV file by cell-type and chromosome-specific. Each row is a sample and each column is a gene.

analysis-runner --access-level test --dataset bioheart --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 --description "pseudobulk" --output-dir "str/associatr/input_files" pseudobulk.py \
--input-dir=gs://cpg-bioheart-test/str/anndata/saige-qtl/anndata_objects_from_HPC --cell-types=CD4_TCM --chromosomes=1 --job-memory=highmem --job-cpu=16

"""
import csv
import logging
import math

import click
import numpy as np
import pandas as pd
import scanpy as sc

from cpg_utils import to_path
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
@click.option('--input-dir', help='GCS Path to the input AnnData object')
@click.option('--sample-id-file-path', help='GCS Path to the sample ID file to do the analysis on')
@click.option('--cell-types', help='Name of the cell type, comma separated if multiple')
@click.option('--chromosomes', help='Chromosome number eg 1, comma separated if multiple')
@click.option('--job-storage', help='Storage of the batch job eg 30G', default='8G')
@click.option('--job-memory', help='Memory of the batch job', default='standard')
@click.option('--job-cpu', help='Number of CPUs of Hail batch job', default=8)
@click.option(
    '--target-sum',
    help='Target sum of counts per cell (for normalization purposes)',
    default=1e6,
)
@click.option(
    '--min-pct',
    help='Minimum percentage of cells expressing a gene to be included in the pseudobulk',
    default=1,
)
@click.command()
def main(
    input_dir,
    sample_id_file_path,
    cell_types,
    chromosomes,
    job_storage,
    job_memory,
    job_cpu,
    target_sum,
    min_pct,
):
    """
    Perform pseudobulk (mean aggregation) of an input AnnData object

    """
    b = get_batch('Run pseudobulk')

    logging.info(f'Cell types to run: {cell_types}')
    logging.info(f'Chromosomes to run: {chromosomes}')

    for cell_type in cell_types.split(','):
        for chromosome in chromosomes.split(','):
            input_file = f'{input_dir}/{cell_type}_chr{chromosome}.h5ad'
            j = b.new_python_job(name=f'Pseudobulk for {cell_type}: chr{chromosome}')
            j.image(image_path('scanpy'))
            j.cpu(job_cpu)
            j.memory(job_memory)
            j.storage(job_storage)
            j.call(pseudobulk, input_file, sample_id_file_path, target_sum, min_pct)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
