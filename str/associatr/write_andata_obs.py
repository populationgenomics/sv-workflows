#!/usr/bin/env python3

"""
This script writes the andata obs to a gcs bucket

analysis-runner --dataset "bioheart" --access-level "test" --description "get cell counts per individual" --output-dir "str/associatr/helper" --image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
--storage=20G --memory=standard --cpu=16 \
write_andata_obs.py --input-h5ad-file-path=gs://cpg-bioheart-test/str/anndata/concatenated_gene_info_donor_info_obs.csv
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
    #create a dataframe with CPG ID, cell type, and counts of cells
    counts = df.groupby(['cpg_id', 'wg2_scpred_prediction']).size()

    input_file_name = input_h5ad_file_path.split('/')[-1].split('.')[0]

    output_gcs = output_path(
        f'{version}/{input_file_name}_cell_counts_per_individual.csv', 'analysis'
    )

    # write to CSV
    counts.to_csv(output_gcs, header=['count'])


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
