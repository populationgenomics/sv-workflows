#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments
"""
This script runs associaTR, given a chromosome or cell type.
Ensure prior scripts have been run to generate dependent files, particularly:
- get_cis_numpy_files.py
- pseudobulk.py
- get_covariates.py
- qc_filters_associatr.py (depends on qc_annotator.py)

 analysis-runner --dataset "bioheart" --config associatr_runner.toml \
    --description "run associatr" \
    --access-level "test" \
    --output-dir "str/associatr/tob_n1055" \
     python3 associatr_runner.py


"""
import json

import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def gene_cis_window_file_reader(file_path):
    cis_window = pd.read_csv(file_path, sep='\t', header=None)
    cis_window.columns = ['chrom', 'start', 'end']
    return f'{cis_window["chrom"][0]}:{cis_window["start"][0]}-{cis_window["end"][0]}'


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

    for celltype in get_config()['associatr']['celltypes'].split(','):
        for chromosome in get_config()['associatr']['chromosomes'].split(','):
            input_dir = get_config()['associatr']['vcf_file_dir']
            vcf_file_path = f'{input_dir}/hail_filtered_{chromosome}.vcf.bgz'
            variant_vcf = b.read_input_group(
                **dict(
                    base=vcf_file_path,
                    tbi=vcf_file_path + '.tbi',
                ),
            )
            gene_list_dir = get_config()['associatr']['gene_list_dir']
            with to_path(f'{gene_list_dir}/{celltype}/{chromosome}_{celltype}_gene_list.json').open('r') as file:
                pseudobulk_gene_names = json.load(file)
            for gene in pseudobulk_gene_names:
                cis_window_dir = get_config()['associatr']['cis_window_dir']
                cis_window_size = get_config()['associatr']['cis_window_size']
                version = get_config()['associatr']['version']

                if to_path(
                    output_path(
                        f'results/{version}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.tsv',
                        'analysis',
                    ),
                ).exists():
                    continue

                gene_cis_window_file = f'{cis_window_dir}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.bed'
                # need to extract the gene start and end from the cis window file for input into 'region'
                cis_window_region = gene_cis_window_file_reader(gene_cis_window_file)
                pheno_cov_numpy_dir = get_config()['associatr']['pheno_cov_numpy_dir']
                gene_pheno_cov = b.read_input(f'{pheno_cov_numpy_dir}/{celltype}/{chromosome}/{gene}_pheno_cov.npy')

                # run associaTR job on the gene
                associatr_job = b.new_job(name=f'Run associatr on {gene} [{celltype};{chromosome}]')
                if get_config()['associatr']['always_run']:
                    associatr_job.always_run()

                associatr_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/trtools_nonlinear_sum:piecewise-reg')
                associatr_job.storage(get_config()['associatr']['job_storage'])
                associatr_job.cpu(get_config()['associatr']['job_cpu'])
                associatr_job.declare_resource_group(association_results={'tsv': '{root}.tsv'})
                associatr_job.command(
                    f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.base} {celltype}_{chromosome}_{gene} {gene_pheno_cov} --region={cis_window_region} --vcftype=eh",
                )
                b.write_output(
                    associatr_job.association_results,
                    output_path(
                        f'results/{version}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp',
                        'analysis',
                    ),
                )
                manage_concurrency_for_job(associatr_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
