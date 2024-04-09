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
    --output-dir "str/associatr/240_libraries_tenk10kp1_v2_run/genome-wide-run/cohort_1" \
     python3 associatr_runner.py


"""
import json
import pandas as pd
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path, init_batch




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
            job.depends_on(
                _dependent_jobs[-get_config()['associatr']['max_parallel_jobs']]
            )
        _dependent_jobs.append(job)

    for celltype in get_config()['associatr']['celltypes'].split(','):
        vcf_file_path = f'gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1/hail_filtered.vcf.bgz'
        variant_vcf = b.read_input_group(
            **dict(
                base=vcf_file_path,
                tbi=vcf_file_path + '.tbi',
            )
        )

        for gene in ['ENSG00000269893','ENSG00000111666','ENSG00000164308','ENSG00000135968']:
            if gene == 'ENSG00000269893':
                chromosome = 'chr4'
            elif gene == 'ENSG00000111666':
                chromosome = 'chr12'
            elif gene == 'ENSG00000164308':
                chromosome = 'chr5'
            else:
                chromosome = 'chr2'
            pheno_cov_numpy_dir = get_config()['associatr']['pheno_cov_numpy_dir']
            gene_pheno_cov = b.read_input(
                f'{pheno_cov_numpy_dir}/{celltype}/{chromosome}/{gene}_pheno_cov.npy'
            )

            # run associaTR job on the gene
            associatr_job = b.new_job(
                name=f'Run associatr on {gene} [{celltype};{chromosome}]'
            )
            associatr_job.image(get_config()['images']['trtools'])
            associatr_job.storage(get_config()['associatr']['job_storage'])
            associatr_job.cpu(get_config()['associatr']['job_cpu'])
            associatr_job.declare_resource_group(
                association_results={'tsv': '{root}.tsv'}
            )
            associatr_job.command(
                f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.base} {celltype}_{chromosome}_{gene} {gene_pheno_cov} --vcftype=eh"
            )
            version = get_config()['associatr']['version']
            b.write_output(
                associatr_job.association_results,
                output_path(
                    f'results/{version}/{celltype}/{chromosome}/{gene}_genome_wide_bp',
                    'analysis',
                ),
            )
            manage_concurrency_for_job(associatr_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
