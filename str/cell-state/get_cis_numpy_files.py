#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script aims to:
 - output gene lists for each cell type, chromosome, pathway combination (after filtering out lowly expressed genes)
 - output cis window files for each gene with scRNA data (cell type + chr specific)
 - optionally removes samples based on a provided sample file
 - perform rank-based inverse normal transformation on pseudobulk data (per gene basis)
 - output gene-level phenotype and covariate numpy objects for input into eqtl pipeline

 analysis-runner  --config get_cis_numpy_files.toml --dataset "tenk10k" --access-level "test" \
--description "get cis and numpy" --output-dir "str/cellstate/input_files/tob" \
--image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
python3 get_cis_numpy_files.py

"""
import json
from ast import literal_eval

import numpy as np
import pandas as pd
import scanpy as sc
from cyvcf2 import VCF
from scipy.stats import norm

import hail as hl
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, image_path, init_batch, output_path



def cis_window_numpy_extractor(
    input_h5ad_dir,
    input_pseudobulk_dir,
    input_cov_dir,
    chromosome,
    cell_type,
    cis_window,
    version,
    chrom_len,
    min_pct,
    pathway
):
    """
    Creates gene-specific cis window files and phenotype-covariate numpy objects

    """
    init_batch()

    # read in anndata object because anndata.vars has the start, end coordinates of each gene
    h5ad_file_path = f'{input_h5ad_dir}/{cell_type}_{chromosome}.h5ad'
    adata = sc.read_h5ad(to_path(h5ad_file_path))

    # read in covariate files
    covariate_path = f'{input_cov_dir}/{cell_type}_covariates.csv'
    covariates = pd.read_csv(covariate_path)

    master_pseudobulk=[]
    #read in pseudobulk data
    pseudobulk_low = pd.read_csv(f'{input_pseudobulk_dir}/{cell_type}/{cell_type}_{chromosome}_low_{pathway}_pseudobulk.csv')
    pseudobulk_medium = pd.read_csv(f'{input_pseudobulk_dir}/{cell_type}/{cell_type}_{chromosome}_medium_{pathway}_pseudobulk.csv')
    pseudobulk_high = pd.read_csv(f'{input_pseudobulk_dir}/{cell_type}/{cell_type}_{chromosome}_high_{pathway}_pseudobulk.csv')
    pseudobulk_low['activity'] = 'low'  # add activity column to pseudobulk df
    pseudobulk_medium['activity'] = 'medium'  # add activity column to pseudobulk df
    pseudobulk_high['activity'] = 'high'  # add activity

    #find the common columns (ie genes) in all pseudobulk dataframes
    # Find common columns in the order they appear in df1
    common_cols = [col for col in pseudobulk_low.columns if col in pseudobulk_medium.columns and col in pseudobulk_high.columns]

    # Subset all DataFrames using this ordered list
    df1_common = pseudobulk_low[common_cols]
    df2_common = pseudobulk_medium[common_cols]
    df3_common = pseudobulk_high[common_cols]

    # Concatenate the dataframes along the rows
    master_pseudobulk.append(df1_common)
    master_pseudobulk.append(df2_common)
    master_pseudobulk.append(df3_common)
    # Concatenate all activities into a single DataFrame
    pseudobulk = pd.concat(master_pseudobulk, axis=0)

    # extract genes in pseudobulk df
    gene_names = [col for col in pseudobulk.columns if col not in ["individual", "activity"]]

    # write filtered gene names to a JSON file
    with to_path(
        output_path(
            f'scRNA_gene_lists/{pathway}/{min_pct}_min_pct_cells_expressed/{cell_type}/{chromosome}_{cell_type}_gene_list.json',
        ),
    ).open('w') as write_handle:
        json.dump(gene_names, write_handle)

    for gene in gene_names:
        # get gene body position (start and end) and add window
        start_coord = adata.var[adata.var.index == gene]['start']
        end_coord = adata.var[adata.var.index == gene]['end']

        left_boundary = max(1, int(start_coord.iloc[0]) - int(cis_window))
        right_boundary = min(int(end_coord.iloc[0]) + int(cis_window), chrom_len)

        data = {'chromosome': chromosome, 'start': left_boundary, 'end': right_boundary}
        ofile_path = output_path(f'cis_window_files/{pathway}/{cell_type}/{chromosome}/{gene}_{cis_window}bp.bed')
        # write cis window file to gcp
        pd.DataFrame(data, index=[gene]).to_csv(ofile_path, sep='\t', header=False, index=False)

        # make the phenotype-covariate numpy objects
        pseudobulk.rename(columns={'individual': 'sample_id'}, inplace=True)  # noqa: PD002
        gene_pheno = pseudobulk[['sample_id', 'activity',gene]]

        # rank-based inverse normal transformation based on R's orderNorm()
        # Rank the values
        gene_pheno.loc[:, 'gene_rank'] = gene_pheno[gene].rank()
        # Calculate the percentile of each rank
        gene_pheno.loc[:, 'gene_percentile'] = (gene_pheno.loc[:, 'gene_rank'] - 0.5) / (len(gene_pheno))
        # Use the inverse normal cumulative distribution function (quantile function) to transform percentiles to normal distribution values
        gene_pheno.loc[:, 'gene_inverse_normal'] = norm.ppf(gene_pheno.loc[:, 'gene_percentile'])
        gene_pheno = gene_pheno[['sample_id', 'activity', 'gene_inverse_normal']]


        gene_pheno_cov = gene_pheno.merge(covariates, on='sample_id')

        # filter for samples that were assigned a CPG ID; unassigned samples after demultiplexing will not have a CPG ID
        gene_pheno_cov = gene_pheno_cov[gene_pheno_cov['sample_id'].str.startswith('CPG')]

        # write out as a csv
        gene_pheno_cov.to_csv(
            output_path(f'pheno_cov_csv/{pathway}/{cell_type}/{chromosome}/{gene}_pheno_cov.csv'),
            index=False,
        )


def main():
    """
    Run cis window extraction and phenotype/covariate numpy object creation
    """
    b = get_batch(name='get cis_numpy files')

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= get_config()['get_cis_numpy']['max_parallel_jobs']:
            job.depends_on(_dependent_jobs[-get_config()['get_cis_numpy']['max_parallel_jobs']])
        _dependent_jobs.append(job)

    for cell_type in get_config()['get_cis_numpy']['cell_types'].split(','):
        for chrom in get_config()['get_cis_numpy']['chromosomes'].split(','):
            init_batch()
            chrom_len = hl.get_reference('GRCh38').lengths[chrom]
            j = b.new_python_job(
                name=f'Extract cis window & phenotype and covariate numpy object for {cell_type}: {chrom}',
            )
            j.image(image_path('scanpy'))
            j.cpu(get_config()['get_cis_numpy']['job_cpu'])
            j.memory(get_config()['get_cis_numpy']['job_memory'])
            j.storage(get_config()['get_cis_numpy']['job_storage'])

            j.call(
                cis_window_numpy_extractor,
                get_config()['get_cis_numpy']['input_h5ad_dir'],
                get_config()['get_cis_numpy']['input_pseudobulk_dir'],
                get_config()['get_cis_numpy']['input_cov_dir'],
                chrom,
                cell_type,
                get_config()['get_cis_numpy']['cis_window'],
                get_config()['get_cis_numpy']['version'],
                chrom_len,
                get_config()['get_cis_numpy']['min_pct'],
                get_config()['get_cis_numpy']['pathway'],
            )

            manage_concurrency_for_job(j)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments
