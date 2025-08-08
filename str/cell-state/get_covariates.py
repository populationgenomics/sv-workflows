#!/usr/bin/env python3
"""
This script calculates cell-type specific PCs from the pseudobulk RNA data (genome-wide),
merges them with other pre-calculated covariates, and writes the file to GCP.


analysis-runner --access-level test --dataset tenk10k --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
--description "Get covariates" --output-dir "str/cellstate/input_files/meanpool/stratified/tob" get_covariates.py --input-dir=gs://cpg-tenk10k-test/str/cellstate/input_files/meanpool/stratified/tob/pseudobulk \
--cell-types=B_intermediate  --covariate-file-path=gs://cpg-tenk10k-test/str/associatr/final-freeze/input_files/tob_n950_sample_covariates.csv \
--num-pcs=6 --pathway=GOBP_STEROID_HORMONE_MEDIATED_SIGNALING_PATHWAY_subtype

"""

import logging

import click
import pandas as pd
import scanpy as sc

import hail as hl
import re

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def get_covariates(pseudobulk_input_dir, cell_type, covariate_file_path, num_pcs, pathway,activity_level):
    """
    Calculates cell-type specific PCs from the pseudobulk data (genome-wide),
    merges them with other pre-calculated covariates, and writes file to GCP
    """
    init_batch()

    # read in and concatenate pseudobulk data (chr-specific --> genome-wide anndata)
    adatas = []
    gcs_files = list(to_path(f'{pseudobulk_input_dir}/{cell_type}').rglob(f'*{pathway}_pseudobulk.csv'))

    for file in gcs_files:
        # Extract activity label from filename
        match = re.search(r"(high|medium|low)", str(file), re.IGNORECASE)
        activity = match.group(1).lower() if match else "unknown"
        if activity == activity_level:
            adata = sc.read_csv(file)
            adatas.append(adata)
    adata_genome = sc.concat(adatas, axis=1)

    # subset to high variance genes prior to PCA
    sc.pp.highly_variable_genes(adata_genome, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_genome = adata_genome[:, adata_genome.var.highly_variable]

    # unit variance scaling for PCA
    sc.pp.scale(adata_genome, max_value=10)

    # PCA
    sc.tl.pca(adata_genome, svd_solver='arpack')

    # Variance ratio plot
    sc.pl.pca_variance_ratio(adata_genome, save='local.png')
    hl.hadoop_copy(
        'figures/pca_variance_ratiolocal.png',
        output_path(f'pseudobulk_RNA_PCs_scree_plots/{cell_type}_scree_plot.png'),
    )

    # extract PCs
    df_pcs = pd.DataFrame(adata_genome.obsm['X_pca'])
    df_pcs.index = adata_genome.obs.index
    # index (CPG ids) are not stored in 'sample_id' column
    df_pcs = df_pcs.rename_axis('sample_id').reset_index()
    df_pcs = df_pcs[['sample_id'] + list(range(num_pcs))]
    # rename PC columns: rna_PC{num}
    df_pcs = df_pcs.rename(columns={i: f'rna_PC{i+1}' for i in range(num_pcs)})

    # read in covariates
    cov = pd.read_csv(covariate_file_path)

    merged_df = cov.merge(df_pcs, on='sample_id')

    # write to GCP
    merged_df.to_csv(
        output_path(f'covariates/{num_pcs}_rna_pcs/{pathway}/{activity_level}/{cell_type}_covariates.csv'),
        index=False,
    )


@click.option('--input-dir', help='GCS Path to the input dir storing pseudobulk CSV files')
@click.option('--cell-types', help='Name of the cell type, comma separated if multiple')
@click.option('--job-storage', help='Storage of the batch job eg 30G', default='8G')
@click.option('--job-memory', help='Memory of the batch job', default='standard')
@click.option('--job-cpu', help='Number of CPUs of Hail batch job', default=8)
@click.option('--covariate-file-path', help='GCS Path to the existing covariate file')
@click.option('--num-pcs', help='Number of RNA PCs to calculate', default=20)
@click.option(
    '--pathway', help='GOBP Pathway name for the analysis', default='GOBP_MULTI_MULTICELLULAR_ORGANISM_PROCESS_subtype'
)
@click.command()
def main(input_dir, cell_types, job_storage, job_memory, job_cpu, covariate_file_path, num_pcs, pathway):
    """
    Obtain cell-type specific covariates for pseudobulk associaTR model
    """
    b = get_batch(name=f'Get covariates for {cell_types}')

    logging.info(f'Cell types to run: {cell_types}')

    for cell_type in cell_types.split(','):
        for activity_level in ['low', 'medium', 'high']:
            # call the get_covariates function
            j = b.new_python_job(name=f'Obtain covariates for for {cell_type} {activity_level}')
            j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3')
            j.cpu(job_cpu)
            j.memory(job_memory)
            j.storage(job_storage)
            j.call(get_covariates, input_dir, cell_type, covariate_file_path, num_pcs, pathway,activity_level)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
