#!/usr/bin/env python3

"""
This script writes the andata obs to a gcs bucket

analysis-runner --dataset "bioheart" --access-level "test" --description "get cell counts per individual" --output-dir "str/associatr/helper" --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
--storage=5G --memory=standard --cpu=4 \
get_cpg_individuals.py --input-h5ad-file-path=gs://cpg-bioheart-test/str/anndata-240/concatenated_gene_info_donor_info_obs.csv
"""

import click
import scanpy as sc
from cpg_utils.hail_batch import output_path
import pandas as pd
from cpg_utils import to_path


@click.option('--input-h5ad-file-path', help='GCS path to the input anndata object')
@click.option('--version', help='version of the output', default='v1')
@click.command()
def main(input_h5ad_file_path: str, version: str):
    """
    Extracts the number of cells per individual per cell type from the anndata object.
    """
    df = pd.read_csv(input_h5ad_file_path)

    length = len(df[df['individual'].str.startswith('CPG')]['individual'].unique())
    length_all = len(df['individual'].unique())

    print(f'Number of individuals with CPG ids:{length}')
    print(f'Number of individuals with all ids:{length_all}')


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter