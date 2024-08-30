#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments
"""
This script runs associaTR on methylation data.

 analysis-runner --dataset "bioheart" --config associatr_runner.toml \
    --description "run associatr" \
    --access-level "test"
    --output-dir "str/associatr-methylation/bioheart_n25" \
     python3 associatr_runner.py


"""
import json

import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path




def main():
    """
    Run associaTR processing pipeline
    """
    b = get_batch(name='Run associatr')
    init_batch()

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= get_config()['associatr']['max_parallel_jobs']:
            job.depends_on(_dependent_jobs[-get_config()['associatr']['max_parallel_jobs']])
        _dependent_jobs.append(job)


    for chromosome in get_config()['associatr']['chromosomes'].split(','):
        vcf_file_path = 'gs://cpg-bioheart-test/str/associatr/gwas-cell-prop/input_files/fm_estr.vcf'
        variant_vcf = b.read_input_group(
            **dict(
                vcf=vcf_file_path,
            ),
        )
        methyl_dir = get_config()['associatr']['pheno_cov_numpy_dir']
        for site_numpy in list(to_path(methyl_dir).glob(f'*.npy')):
            cis_window_size = 50000
            site_coord = str(site_numpy).split('/')[-1].split('_')[1]
            site = f'{chromosome}_{site_coord}'

            if to_path(
                output_path(
                    f'results/{chromosome}/{site}_{cis_window_size}bp.tsv',
                    'analysis',
                ),
            ).exists():
                continue

            # need to extract the gene start and end from the cis window file for input into 'region'
            cis_window_start = max(0, int(site_coord) - cis_window_size)
            cis_window_end = int(site_coord) + cis_window_size
            cis_window_region = f'{chromosome}:{cis_window_start}-{cis_window_end}'

            pheno_cov_numpy_dir = get_config()['associatr']['pheno_cov_numpy_dir']
            gene_pheno_cov = b.read_input(f'{pheno_cov_numpy_dir}/{chromosome}/{site}_pheno_cov.npy')

            # run associaTR job on the gene
            associatr_job = b.new_job(name=f'Run associatr on {site}')

            associatr_job.image(get_config()['images']['trtools_hope_version'])
            associatr_job.storage(get_config()['associatr']['job_storage'])
            associatr_job.cpu(get_config()['associatr']['job_cpu'])
            associatr_job.declare_resource_group(association_results={'tsv': '{root}.tsv'})
            associatr_job.command(
                f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.vcf} {site} {gene_pheno_cov} --region={cis_window_region} --vcftype=eh",
            )
            b.write_output(
                associatr_job.association_results,
                output_path(
                    f'results/{chromosome}/{site}_{cis_window_size}bp',
                    'analysis',
                ),
            )
            manage_concurrency_for_job(associatr_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
