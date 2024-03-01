#!/usr/bin/env python3

"""
This script aims to:
 - output cis window files for each gene with scRNA data (cell type + chr specific)
 - output gene-level phenotype and covariate numpy objects for input into associatr

 analysis-runner --dataset "bioheart" --access-level "test" --description "get cis and numpy" --output-dir "str/associatr/input_files" --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
 get_cis_numpy_files.py --input-h5ad-dir=gs://cpg-bioheart-test/str/anndata/saige-qtl/anndata_objects_from_HPC --input-pseudobulk-dir=gs://cpg-bioheart-test/str/associatr/input_files/pseudobulk/pseudobulk --input-cov-dir=gs://cpg-bioheart-test/str/associatr/input_files/covariates \
 --chromosomes=chr1 --cell-types=ASDC --cis-window=100000 --version=v1

"""

import click
import numpy as np
import pandas as pd
import scanpy as sc
import hail as hl

from cpg_utils.hail_batch import get_batch, output_path
from cpg_utils import to_path


def cis_window_numpy_extractor(
    input_h5ad_dir,
    input_pseudobulk_dir,
    input_cov_dir,
    chromosome,
    cell_type,
    cis_window,
    version,
):
    """
    Creates gene-specific cis window files and phenotype-covaraite numpy objects

    """
    # read in anndata object because anndata.vars has the start, end coordinates of each gene
    h5ad_file_path = f'{input_h5ad_dir}/{cell_type}_{chromosome}.h5ad'
    expression_h5ad_path = to_path(h5ad_file_path).copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path)

    # read in pseudobulk and covariate files
    pseudobulk_path = (
        f'{input_pseudobulk_dir}/{cell_type}/{cell_type}_{chromosome}_pseudobulk.csv'
    )
    pseudobulk = pd.read_csv(pseudobulk_path)
    covariate_path = f'{input_cov_dir}/{cell_type}_covariates.csv'
    covariates = pd.read_csv(covariate_path)

    for gene in adata.var.index:
        # get gene body position (start and end) and add window
        start_coord = adata.var[adata.var.index == gene]['start']
        end_coord = adata.var[adata.var.index == gene]['end']

        left_boundary = max(1, int(start_coord) - cis_window)
        right_boundary = min(
            int(end_coord) + cis_window, hl.get_reference('GRCh38').lengths[chromosome]
        )

        data = {'chromosome': chromosome, 'start': left_boundary, 'end': right_boundary}
        ofile_path = output_path(
            f'cis_window_files/{version}/{cell_type}/{chromosome}/{gene}_{cis_window}bp.bed'
        )
        # write cis window file to gcp
        pd.DataFrame(data, index=[gene]).to_csv(
            ofile_path, sep='\t', header=False, index=False
        )

        # make the phenotype-covariate numpy objects
        pseudobulk.rename(columns={'individual': 'sample_id'}, inplace=True)
        gene_pheno = pseudobulk[['sample_id', gene]]
        gene_pheno_cov = gene_pheno.merge(covariates, on='sample_id', how='inner')
        gene_pheno_cov['sample_id'] = gene_pheno_cov['sample_id'].str[
            3:
        ]  # remove CPG prefix because associatr expects id to be numeric
        gene_pheno_cov = gene_pheno_cov.to_numpy()
        np.save(
            output_path(
                f'pheno_cov_numpy/{version}/{cell_type}/{chromosome}/{gene}_pheno_cov.npy'
            ),
            gene_pheno_cov,
        )


@click.option('--input-h5ad-dir', help='GCS directory to the input AnnData objects')
@click.option(
    '--input-pseudobulk-dir', help='GCS directory to the input pseudobulk CSV files'
)
@click.option('--input-cov-dir', help='GCS directory to the input covariate CSV files')
@click.option('--chromosomes', help=' eg chr22, comma separated if multiple')
@click.option('--cell-types', help='cell type, comma separated if multiple')
@click.option('--cis-window', help='cis window size in bp')
@click.option('--version', default='v1', help='version of the output files')
@click.option('--job-cpu', default=8, help='Number of CPUs to use for the job')
@click.option('--job-memory', default='standard', help='Memory to use for the job')
@click.option('--job-storage', default='10G', help='Storage to use for the job')
@click.command()
def main(
    input_h5ad_dir,
    input_pseudobulk_dir,
    input_cov_dir,
    chromosomes,
    cell_types,
    cis_window,
    version,
    job_cpu,
    job_memory,
    job_storage,
):
    """
    Run cis window extraction and phenotype/covariate numpy object creation
    """
    b = get_batch()
    for cell_type in cell_types.split(','):
        for chrom in chromosomes.split(','):
            j = b.new_python_job(
                name=f'Extract cis window & phenotype and covariate numpy object for {cell_type}: {chrom}'
            )
            j.image(
                'australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3'
            )
            j.cpu(job_cpu)
            j.memory(job_memory)
            j.storage(job_storage)
            j.call(
                cis_window_numpy_extractor,
                input_h5ad_dir,
                input_pseudobulk_dir,
                input_cov_dir,
                chrom,
                cell_type,
                cis_window,
                version,
            )
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
