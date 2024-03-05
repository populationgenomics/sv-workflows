#!/usr/bin/env python3

"""
This script aims to extract the number of cells per individual per cell type from the anndata object (containing all individuals, all chromosomes).
The output is a CSV file with three columns: CPG ID, cell type, counts of cells.

analysis-runner --dataset "bioheart" --access-level "test" --description "get cell counts per individual" --output-dir "str/associatr/helper" --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
--storage=10G --cpu=4 \
cell_per_individual_extractor.py --input-h5ad-dir=gs://cpg-bioheart-test/str/anndata/B_naive_chr10.h5ad
"""

import click
import scanpy as sc
from cpg_utils.hail_batch import output_path
from cpg_utils import to_path


@click.option('--input-h5ad-dir', help='GCS path to the input anndata object', type=str)
@click.option('--version', help='version of the output', default='v1')
@click.command()
def main(input_h5ad_dir, version):
    """
    Extracts the number of cells per individual per cell type from the anndata object.
    """
    # read in anndata object because anndata.obs has the CPG ID and cell type
    expression_h5ad_path = to_path(input_h5ad_dir).copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path)

    # create a dataframe with CPG ID, cell type, and counts of cells
    counts = adata.obs.groupby(['cpg_id', 'wg2_scpred_prediction']).size()

    input_file_name = input_h5ad_dir.split('/')[-1].split('.')[0]

    output_gcs = output_path(
        f'{version}/{input_file_name}_cell_counts_per_individual.csv', 'analysis'
    )

    # write to CSV
    counts.to_csv(output_gcs, header=['count'])


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
