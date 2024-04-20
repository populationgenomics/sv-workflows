#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script aims to:
 - output gene lists for each cell type and chromosome (after filtering out lowly expressed genes)
 - output cis window files for each gene with scRNA data (cell type + chr specific)
 - optionally removes samples based on a provided sample file
 - perform rank-based inverse normal transformation on pseudobulk data (per gene basis)
 - output gene-level phenotype and covariate numpy objects for input into associatr

 analysis-runner  --config get_cis_numpy_files.toml --dataset "bioheart" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr/rna_pc_calibration/2_pcs/input_files"  \
--image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
python3 get_cis_numpy_files.py

"""
import json
from ast import literal_eval

import numpy as np
import pandas as pd
import scanpy as sc
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
    remove_samples_file,
):
    """
    Creates gene-specific cis window files and phenotype-covariate numpy objects

    """
    init_batch()

    # read in anndata object because anndata.vars has the start, end coordinates of each gene
    h5ad_file_path = f'{input_h5ad_dir}/{cell_type}_{chromosome}.h5ad'
    expression_h5ad_path = to_path(h5ad_file_path).copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path)

    # read in pseudobulk and covariate files
    pseudobulk_path = f'{input_pseudobulk_dir}/{cell_type}/{cell_type}_{chromosome}_pseudobulk.csv'
    pseudobulk = pd.read_csv(pseudobulk_path)
    covariate_path = f'{input_cov_dir}/{cell_type}_covariates.csv'
    covariates = pd.read_csv(covariate_path)

    # extract genes in pseudobulk df
    gene_names = list(pseudobulk.columns[1:])  # individual ID is the first column

    # write filtered gene names to a JSON file
    with to_path(
        output_path(
            f'scRNA_gene_lists/{min_pct}_min_pct_cells_expressed/{cell_type}/{chromosome}_{cell_type}_gene_list.json',
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
        ofile_path = output_path(f'cis_window_files/{version}/{cell_type}/{chromosome}/{gene}_{cis_window}bp.bed')
        # write cis window file to gcp
        pd.DataFrame(data, index=[gene]).to_csv(ofile_path, sep='\t', header=False, index=False)

        # make the phenotype-covariate numpy objects
        pseudobulk.rename(columns={'individual': 'sample_id'}, inplace=True)  # noqa: PD002
        gene_pheno = pseudobulk[['sample_id', gene]]

        # remove samples that are in the remove_samples_file
        if remove_samples_file:
            with to_path(remove_samples_file).open() as f:
                array_string = f.read().strip()
                remove_samples = literal_eval(array_string)
                gene_pheno = gene_pheno[~gene_pheno['sample_id'].isin(remove_samples)]

        # rank-based inverse normal transformation based on R's orderNorm()
        # Rank the values
        gene_pheno.loc[:, 'gene_rank'] = gene_pheno[gene].rank()
        # Calculate the percentile of each rank
        gene_pheno.loc[:, 'gene_percentile'] = (gene_pheno.loc[:, 'gene_rank'] - 0.5) / (len(gene_pheno))
        # Use the inverse normal cumulative distribution function (quantile function) to transform percentiles to normal distribution values
        gene_pheno.loc[:, 'gene_inverse_normal'] = norm.ppf(gene_pheno.loc[:, 'gene_percentile'])
        gene_pheno = gene_pheno[['sample_id', 'gene_inverse_normal']]

        gene_pheno_cov = gene_pheno.merge(covariates, on='sample_id', how='inner')

        # filter for samples that were assigned a CPG ID; unassigned samples after demultiplexing will not have a CPG ID
        gene_pheno_cov = gene_pheno_cov[gene_pheno_cov['sample_id'].str.startswith('CPG')]


        gene_pheno_cov['sample_id'] = gene_pheno_cov['sample_id'].str[
            3:
        ]  # remove CPG prefix because associatr expects id to be numeric

        gene_pheno_cov['sample_id'] = gene_pheno_cov['sample_id'].astype(float)

        gene_pheno_cov = gene_pheno_cov.to_numpy()
        with hl.hadoop_open(
            output_path(f'pheno_cov_numpy/{version}/{cell_type}/{chromosome}/{gene}_pheno_cov.npy'),
            'wb',
        ) as f:
            np.save(f, gene_pheno_cov)


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
                get_config()['get_cis_numpy']['remove_samples_file'],
            )

            manage_concurrency_for_job(j)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments
