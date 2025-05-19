#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments
"""
This script runs associaTR on atac data.

 analysis-runner --dataset "bioheart" --config associatr_runner.toml \
    --description "run associatr" \
    --access-level "test" \
    --output-dir "str/associatr-atac/tob/input_files/10kb_estrs/v1" \
     python3 associatr_runner.py


"""
import json

import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path, reset_batch


def main():
    """
    Run associaTR processing pipeline
    """

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

    for cell_type in get_config()['associatr']['cell_types'].split(','):
        vcf_file_dir = get_config()['associatr']['vcf_dir']
        vcf_file_path = f'{vcf_file_dir}/{cell_type}_estrs.sorted.vcf.gz'
        pheno_dir = get_config()['associatr']['pheno_cov_numpy_dir']
        site_numpy_list = list(to_path(f'{pheno_dir}/{cell_type}/pheno_cov_numpy').glob('*.npy'))
        for i in range(0, len(site_numpy_list), 2000):
            _dependent_jobs = []
            reset_batch()
            batch_gene_files = site_numpy_list[i : i + 2000]
            b = get_batch(name='Run associatr-atac' + cell_type + ' sites ' + str(i) + '-' + str(i + 2000))
            variant_vcf = b.read_input_group(
                **dict(
                    vcf=vcf_file_path,
                    vcf_index=f'{vcf_file_path}.tbi',
                ),
            )
            for site_numpy in batch_gene_files:
                cis_window_size = get_config()['associatr']['cis_window_size']
                site_coord = str(site_numpy).split('/')[-1].split('_')[0]
                start_coord = site_coord.split(':')[1].split('-')[0]
                chromosome = site_coord.split(':')[0]
                end_coord = site_coord.split(':')[1].split('-')[1]
                site = site_coord

                if to_path(
                    output_path(
                        f'results/{cell_type}/{site}_{cis_window_size}bp.tsv',
                        'analysis',
                    ),
                ).exists():
                    continue

                # need to extract the gene start and end from the cis window file for input into 'region'
                cis_window_start = max(0, int(start_coord) - cis_window_size)
                cis_window_end = int(end_coord) + cis_window_size
                cis_window_region = f'{chromosome}:{cis_window_start}-{cis_window_end}'

                pheno_cov_numpy_dir = get_config()['associatr']['pheno_cov_numpy_dir']
                gene_pheno_cov = b.read_input(f'{pheno_cov_numpy_dir}/{cell_type}/pheno_cov_numpy/{site}_pheno_cov.npy')

                # run associaTR job on the gene
                associatr_job = b.new_job(name=f'Run associatr on {site}')

                associatr_job.image(get_config()['images']['trtools_hope_version'])
                associatr_job.storage(get_config()['associatr']['job_storage'])
                associatr_job.cpu(get_config()['associatr']['job_cpu'])
                associatr_job.declare_resource_group(association_results={'tsv': '{root}.tsv'})
                associatr_job.command(
                    f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.vcf} {site} {gene_pheno_cov} --region={cis_window_region} --vcftype=eh --non-major-cutoff=0",
                )
                b.write_output(
                    associatr_job.association_results,
                    output_path(
                        f'results/{cell_type}/{site}_{cis_window_size}bp',
                        'analysis',
                    ),
                )
                manage_concurrency_for_job(associatr_job)
            b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter