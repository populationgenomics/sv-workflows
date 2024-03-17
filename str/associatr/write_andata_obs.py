#!/usr/bin/env python3

"""
This script writes the andata obs to a gcs bucket

analysis-runner --dataset "bioheart" --access-level "test" --description "get cell counts per individual" --output-dir "str/associatr/helper/anndata-240" --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
--storage=200G --memory=highmem --cpu=16 \
write_andata_obs.py --input-h5ad-dir=gs://cpg-bioheart-test/str/anndata-240/240_libraries_concatenated_gene_info.h5ad
"""

import click
import scanpy as sc
from cpg_utils.hail_batch import output_path
from cpg_utils import to_path


@click.option('--input-h5ad-file-path', help='GCS path to the input anndata object')
@click.option('--version', help='version of the output', default='v1')
@click.command()
def main(input_h5ad_file_path: str, version: str):
    """
    Extracts the number of cells per individual per cell type from the anndata object.
    """
    # read in anndata object because anndata.obs has the CPG ID and cell type
    expression_h5ad_path = to_path(input_h5ad_file_path).copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path, backed = 'r')

    output_gcs = 'gs://cpg-bioheart-test/str/anndata-240/concatenated_gene_info_donor_info_obs.csv'
    # write to CSV
    adata.obs.to_csv(output_gcs)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter