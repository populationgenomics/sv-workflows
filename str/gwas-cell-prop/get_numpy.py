#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script aims to:
 - optionally removes samples based on a provided sample file
 - perform rank-based inverse normal transformation on cell proportions data
 - output phenotype and covariate numpy objects for input into associatr

 analysis-runner  --dataset "bioheart" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr/gwas-cell-prop/tob_n1055/input_files" \
 --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
 get_numpy.py --covariate-file-path=gs://cpg-bioheart-test/str/associatr/tob_n1055/input_files/tob_covariates_str_run_v1.csv \
    --cell-summary-file-path=gs://cpg-bioheart-test/str/associatr/gwas-cell-prop/input_files/cell_prop.csv \
    --remove-samples-file=gs://cpg-bioheart-test/str/associatr/input_files/remove-samples.txt

"""

import click

from cpg_utils.hail_batch import get_batch


def get_numpy(macrogroup, cov_path, cell_summary_file, remove_samples_file):
    import numpy as np
    from ast import literal_eval
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path
    import hail as hl


    import pandas as pd
    from scipy.stats import norm
    # read in cellsummary file and RINT
    cell_summary = pd.read_csv(cell_summary_file)
    cell_summary = cell_summary[cell_summary['macro_group'] == macrogroup]
    cell_summary = cell_summary[['cpg_id', 'cell_proportion']]

    # remove any samples if specified
    with to_path(remove_samples_file).open() as f:
        array_string = f.read().strip()
        remove_samples = literal_eval(array_string)
        cell_summary = cell_summary[~cell_summary['sample_id'].isin(remove_samples)]

    # RINT
    cell_summary.loc[:, 'rank'] = cell_summary['cell_proportion'].rank()
    # Calculate the percentile of each rank
    cell_summary.loc[:, 'percentile'] = (cell_summary.loc[:, 'rank'] - 0.5) / (len(cell_summary))
    # Use the inverse normal cumulative distribution function (quantile function) to transform percentiles to normal distribution values
    cell_summary.loc[:, 'inverse_normal'] = norm.ppf(cell_summary.loc[:, 'percentile'])
    cell_summary = cell_summary[['cpg_id', 'inverse_normal']]

    # merge with other covariates
    covariates = pd.read_csv(cov_path, sep='\t')
    df = pd.merge(cell_summary, covariates, left_on='cpg_id', right_on='sample_id')
    df = df.drop(columns=['sample_id'])

    # filter for samples with a CPG id
    df = df[df['cpg_id'].str.startswith('CPG')]

    df['cpg_id'] = df['cpg_id'].str[3:]  # remove CPG prefix because associatr expects id to be numeric

    df['cpg_id'] = df['cpg_id'].astype(float)

    df = df.to_numpy()
    with hl.hadoop_open(
        output_path(f'pheno_cov_numpy/{macrogroup}_pheno_cov.npy'),
        'wb',
    ) as f:
        np.save(f, df)


@click.option('--remove-samples-file', required=False, help='Path to the file containing the samples to remove')
@click.option('--covariate-file-path', required=True, help='Path to the covariate file')
@click.option('--cell-summary-file-path', required=True, help='Path to the cell summary file')
@click.command()
def main(covariate_file_path, cell_summary_file_path, remove_samples_file):
    b = get_batch(name='get numpy files')
    for macrogroup in [
        'B',
        'CD4T',
        'CD8T',
        'Cycling',
        'DC',
        'MAIT',
        'Mono',
        'NK',
        'Platelet',
        'Treg',
        'gdT',
        'nan',
        'Eryth',
    ]:
        j = b.new_python_job(name=f'get_numpy_{macrogroup}')
        j.call(get_numpy, macrogroup, covariate_file_path, cell_summary_file_path, remove_samples_file)

    b.run(wait=False)


if __name__ == '__main__':
    main()
