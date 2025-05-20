#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script aims to:
 - output ATAC-seq site-level phenotype and covariate numpy objects for input into associatr per cell type.

 analysis-runner  --config get_cis_numpy.toml --dataset "bioheart" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr-atac/tob/input_files/10kb_estrs" \
python3 get_cis_numpy.py

"""

import pandas as pd

import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch


def cis_window_numpy_extractor(
    batch_sites,
    covariate_path,
    pseudobulk_dir,
    cell_type,
):
    """
    Creates phenotype-covariate numpy objects

    """
    import numpy as np
    import pandas as pd
    from scipy.stats import norm

    import hail as hl

    from cpg_utils import to_path
    from cpg_utils.hail_batch import init_batch, output_path

    init_batch()
    # read in batch_sites and covariate files
    pheno = batch_sites



    for site in pheno['loc']:
        # make the phenotype-covariate numpy objects

        site_pheno = pheno[pheno['loc'] == site]
        site_pheno = site_pheno.drop(columns=['loc'])

        # Now, let's reshape the dataframe
        site_pheno = pd.melt(site_pheno, value_vars=site_pheno.columns, var_name='sample_id', value_name='z_score')

        # match up IDs with CPG IDs
        onek1k_mapper = pd.read_csv('gs://cpg-bioheart-test/str/atac-seq/sample-mapping/OneK1K_sex.tsv', sep='\t')
        cpg_mapper = pd.read_csv('gs://cpg-bioheart-test/str/atac-seq/sample-mapping/saigesaige-qtl_input_files_mapping_for_anna.csv')
        cpg_mapper['external_id'] = cpg_mapper['external_id'].str.replace('-PBMC$', '', regex=True)
        master_mapper = onek1k_mapper.merge(cpg_mapper, left_on = 'ExternalID', right_on = 'external_id')
        df1 = pd.read_csv('gs://cpg-bioheart-test//str/associatr/final-freeze/input_files/tob_n950_sample_covariates.csv')
        mapping_key = df1.merge(master_mapper, left_on = 'sample_id', right_on = 'new_CPG_id')[['sampleid','new_CPG_id']]
        site_pheno = site_pheno.merge(mapping_key, left_on='sample_id', right_on='sampleid', how='inner')
        site_pheno = site_pheno[['new_CPG_id', 'z_score']]
        site_pheno = site_pheno.rename(columns={'new_CPG_id': 'sample_id'})

        # read in RNA covariates
        rna_covariates =pd.read_csv(f'{pseudobulk_dir}/{cell_type}/covariates.txt')
        # Step 1: Filter for rows where 'id' starts with "PC"
        filtered_df = rna_covariates[rna_covariates['id'].str.startswith('PC')]
        transposed_df = filtered_df.set_index('id').T.reset_index()
        transposed_df = transposed_df.rename(columns={'index': 'sample_id'})
        transposed_df = transposed_df.merge(mapping_key, left_on='sample_id', right_on='sampleid', how='inner')
        transposed_df = transposed_df.drop(columns=['sample_id', 'sampleid'])
        print(transposed_df)

        covariates = pd.read_csv(covariate_path)
        # merge with geno PCs, age, sex
        covariates = covariates.merge(transposed_df, left_on = 'sample_id', right_on = 'new_CPG_id').drop(columns=['new_CPG_id'])
        print(covariates)

        print(site_pheno)

        # merge covariates with pseudobulk atac data
        site_pheno_cov = site_pheno.merge(covariates, on='sample_id', how='inner')

        site_pheno_cov['sample_id'] = site_pheno_cov['sample_id'].str[
            3:
        ]  # remove CPG prefix because associatr expects id to be numeric

        site_pheno_cov['sample_id'] = site_pheno_cov['sample_id'].astype(float)

        # shuffle the 'sample_id' column
        site_pheno_cov['sample_id'] = site_pheno_cov['sample_id'].sample(frac=1).reset_index(drop=True)

        site_pheno_cov = site_pheno_cov.to_numpy()
        with hl.hadoop_open(
            output_path(f'{cell_type}_permuted/pheno_cov_numpy/{site}_pheno_cov.npy'),
            'wb',
        ) as f:
            np.save(f, site_pheno_cov)


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

    sites_per_job = 50

    for cell_type in get_config()['get_cis_numpy']['cell_types'].split(','):
        pseudobulk_dir = get_config()['get_cis_numpy']['pseudobulk_dir']
        loc_filtered_dir = get_config()['get_cis_numpy']['loc_filtered_dir']
        pseudobulk = pd.read_csv(f'{pseudobulk_dir}/{cell_type}/PseudobulkMatStandardised.csv')
        #filter down to only the sites in the loc_filtered file (ie sites that are within 10kb of an etr)
        loc_filtered = pd.read_csv(f'{loc_filtered_dir}/{cell_type}_loc.csv')
        subset_df = pseudobulk[pseudobulk['loc'].isin(loc_filtered['loc'])]

        for i in range(0, subset_df.shape[0], sites_per_job):
            batch_sites = subset_df.iloc[i : i + sites_per_job]

            j = b.new_python_job(
                name=f'Extract cis window & phenotype and covariate numpy object: {cell_type}, sites {i}-{i+sites_per_job}',
            )
            j.cpu(get_config()['get_cis_numpy']['job_cpu'])
            j.memory(get_config()['get_cis_numpy']['job_memory'])
            j.storage(get_config()['get_cis_numpy']['job_storage'])

            j.call(
                cis_window_numpy_extractor,
                batch_sites,
                get_config()['get_cis_numpy']['covariate_path'],
                get_config()['get_cis_numpy']['pseudobulk_dir'],
                cell_type,
            )

        manage_concurrency_for_job(j)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments