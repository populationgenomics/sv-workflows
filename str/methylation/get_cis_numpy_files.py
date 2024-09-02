#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script aims to:
 - output list of CpG sites for each chromosome
 - perform rank-based inverse normal transformation on methylation (mod_score) (per CpG site basis)
 - output CpG site-level phenotype and covariate numpy objects for input into associatr

 analysis-runner  --config get_cis_numpy_files.toml --dataset "bioheart" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr-methylation/bioheart_n25/input_files/5kb" \
python3 get_cis_numpy_files.py

"""

import pandas as pd

import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch


def cis_window_numpy_extractor(
    batch_sites,
    covariate_path,
    chromosome,
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
    covariates = pd.read_csv(covariate_path)
    pheno['meth_id'] = pheno['chrom'].astype(str) + '_' + pheno['start'].astype(str)

    for site in pheno['meth_id']:
        # make the phenotype-covariate numpy objects

        site_pheno = pheno[pheno['meth_id'] == site]
        site_pheno = site_pheno.drop(columns=['chrom', 'start', 'meth_id', 'std_dev'])

        # Now, let's reshape the dataframe
        site_pheno = pd.melt(site_pheno, value_vars=site_pheno.columns, var_name='sample_id', value_name=site)

        # rank-based inverse normal transformation based on R's orderNorm()
        # Rank the values
        site_pheno.loc[:, 'gene_rank'] = site_pheno[site].rank()
        # Calculate the percentile of each rank
        site_pheno.loc[:, 'gene_percentile'] = (site_pheno.loc[:, 'gene_rank'] - 0.5) / (len(site_pheno))
        # Use the inverse normal cumulative distribution function (quantile function) to transform percentiles to normal distribution values
        site_pheno.loc[:, 'gene_inverse_normal'] = norm.ppf(site_pheno.loc[:, 'gene_percentile'])
        site_pheno = site_pheno[['sample_id', 'gene_inverse_normal']]

        # match up FS ID with CPG IDs
        mapper = pd.read_csv('gs://cpg-bioheart-test/str/pacbio-methylation-v3/mapper_trgt.csv')
        site_pheno = site_pheno.merge(mapper, left_on='sample_id', right_on='FS_ID', how='inner')
        site_pheno = site_pheno.drop(columns=['sample_id', 'FS_ID'])

        # rename 's' to 'sample_id'
        site_pheno = site_pheno.rename(columns={'s': 'sample_id'})
        site_pheno = site_pheno[['sample_id', 'gene_inverse_normal']]

        site_pheno_cov = site_pheno.merge(covariates, on='sample_id', how='inner')

        site_pheno_cov['sample_id'] = site_pheno_cov['sample_id'].str[
            3:
        ]  # remove CPG prefix because associatr expects id to be numeric

        site_pheno_cov['sample_id'] = site_pheno_cov['sample_id'].astype(float)

        site_pheno_cov = site_pheno_cov.to_numpy()
        with hl.hadoop_open(
            output_path(f'pheno_cov_numpy/{chromosome}/{site}_pheno_cov.npy'),
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

    for chrom in get_config()['get_cis_numpy']['chromosomes'].split(','):
        methyl_dir = get_config()['get_cis_numpy']['input_methyl_dir']
        pheno = pd.read_csv(f'{methyl_dir}/methylation_combined_{chrom}.bed', sep='\t')
        for i in range(0, pheno.shape[0], sites_per_job):
            batch_sites = pheno.iloc[i : i + sites_per_job]

            j = b.new_python_job(
                name=f'Extract cis window & phenotype and covariate numpy object: {chrom}, sites {i}-{i+sites_per_job}',
            )
            j.cpu(get_config()['get_cis_numpy']['job_cpu'])
            j.memory(get_config()['get_cis_numpy']['job_memory'])
            j.storage(get_config()['get_cis_numpy']['job_storage'])

            j.call(
                cis_window_numpy_extractor,
                batch_sites,
                get_config()['get_cis_numpy']['covariate_path'],
                chrom,
            )

        manage_concurrency_for_job(j)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments
