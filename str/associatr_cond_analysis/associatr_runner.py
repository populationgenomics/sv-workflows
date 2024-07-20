#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments
"""
This script runs associaTR (conditioned on the index STR), given a chromosome or cell type.
Ensure prior scripts have been run to generate dependent files, particularly:
- get_cis_numpy_files.py (conditional analysis version)
- pseudobulk.py
- get_covariates.py
- qc_filters_associatr.py (depends on qc_annotator.py)

 analysis-runner --dataset "bioheart" --config associatr_runner_snps.toml \
    --description "run associatr" \
    --access-level "test" \
    --output-dir "str/associatr/cond_analysis/common_variants_snps/tob_n1055" \
     python3 associatr_runner.py

     analysis-runner --dataset "bioheart" --config associatr_runner.toml \
    --description "run associatr" \
    --access-level "test" \
    --output-dir "str/associatr/cond_analysis/bioheart_n990" \
     python3 associatr_runner.py


"""
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

    egenes = pd.read_csv(get_config()['associatr']['egene'])
    celltypes = egenes['cell_type'].unique()

    for cell_type in celltypes:
        egenes_cell = egenes[egenes['cell_type'] == cell_type]
        for chrom in range(1, 23):
            egenes_cell_chrom = egenes_cell[egenes_cell['chr'] == f'chr{chrom}']
            if egenes_cell_chrom.empty:
                continue
            input_dir = get_config()['associatr']['vcf_file_dir']
            vcf_file_path = f'{input_dir}/hail_filtered_chr{chrom}.vcf.bgz'
            variant_vcf = b.read_input_group(
                **dict(
                    base=vcf_file_path,
                    tbi=vcf_file_path + '.tbi',
                ),
            )
            for gene in egenes_cell_chrom['gene_name']:
                cis_window_dir = get_config()['associatr']['cis_window_dir']
                cis_window_size = get_config()['associatr']['cis_window_size']
                version = get_config()['associatr']['version']

                if to_path(
                    output_path(
                        f'results/{version}/{cell_type}/chr{chrom}/{gene}_{cis_window_size}bp.tsv',
                        'analysis',
                    ),
                ).exists():
                    continue

                gene_cis_window_file = f'{cis_window_dir}/{cell_type}/chr{chrom}/{gene}_{cis_window_size}bp.bed'
                # need to extract the gene start and end from the cis window file for input into 'region'
                try:
                    cis_window_region = gene_cis_window_file_reader(gene_cis_window_file)
                except FileNotFoundError:
                    print(f'File not found: {gene_cis_window_file}..skipping')
                    continue
                pheno_cov_numpy_dir = get_config()['associatr']['pheno_cov_numpy_dir']
                gene_pheno_cov = b.read_input(f'{pheno_cov_numpy_dir}/{cell_type}/chr{chrom}/{gene}_pheno_cov.npy')

                # run associaTR job on the gene
                associatr_job = b.new_job(name=f'Run associatr on {gene} [{cell_type};chr{chrom}]')
                if get_config()['associatr']['always_run']:
                    associatr_job.always_run()

                associatr_job.image(get_config()['images']['trtools'])
                associatr_job.storage(get_config()['associatr']['job_storage'])
                associatr_job.cpu(get_config()['associatr']['job_cpu'])
                associatr_job.declare_resource_group(association_results={'tsv': '{root}.tsv'})
                associatr_job.command(
                    f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.base} {cell_type}_chr{chrom}_{gene} {gene_pheno_cov} --region={cis_window_region} --vcftype=eh",
                )
                b.write_output(
                    associatr_job.association_results,
                    output_path(
                        f'results/{version}/{cell_type}/chr{chrom}/{gene}_{cis_window_size}bp',
                        'analysis',
                    ),
                )
                manage_concurrency_for_job(associatr_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
