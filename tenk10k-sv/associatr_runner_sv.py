#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments,too-many-locals
"""
SV-specific associaTR runner.

Difference vs. str/associatr/associatr_runner.py:
- Uses INFO/END from the EH-style SV VCF (produced by
  tenk10k-sv/sv_vcf_for_associatr.py) to include any SV whose interval
  [POS, END] overlaps the gene cis window by >=1 bp, rather than only variants
  whose POS falls inside the window. This is done per gene with a bcftools
  pre-filter before associaTR is invoked:
      bcftools view -r {chrom}:1-{cis_end} -i 'INFO/END>={cis_start}' ...
  associaTR is then given --region {chrom}:1-{chrom_len} so it walks every
  record left in the pre-filtered VCF.
- Assumes the cis window bed files were produced with a 1 Mb window
  (see tenk10k-sv/get_cis_numpy_files_sv.toml).

 analysis-runner --dataset tenk10k --config associatr_runner_sv.toml \
    --description "run associatr on SVs" \
    --access-level test \
    --memory "8G" \
    --output-dir "tenk10k-sv/sv/sensitivity_analysis/tob_n935/1pc" \
    python3 associatr_runner_sv.py
"""

import json

import pandas as pd

import hail as hl
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def gene_cis_window_file_reader(file_path):
    cis_window = pd.read_csv(file_path, sep='\t', header=None)
    cis_window.columns = ['chrom', 'start', 'end']
    return cis_window['chrom'][0], int(cis_window['start'][0]), int(cis_window['end'][0])


def main():
    b = get_batch(name='Run associatr on SVs')
    init_batch()

    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        if len(_dependent_jobs) >= get_config()['associatr']['max_parallel_jobs']:
            job.depends_on(_dependent_jobs[-get_config()['associatr']['max_parallel_jobs']])
        _dependent_jobs.append(job)

    for celltype in get_config()['associatr']['celltypes'].split(','):
        for chromosome in get_config()['associatr']['chromosomes'].split(','):
            chrom_len = hl.get_reference('GRCh38').lengths[chromosome]

            input_dir = get_config()['associatr']['vcf_file_dir']
            vcf_file_path = f'{input_dir}/hail_filtered_{chromosome}.vcf.bgz'
            variant_vcf = b.read_input_group(
                base=vcf_file_path,
                tbi=vcf_file_path + '.tbi',
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
                cis_chrom, cis_start, cis_end = gene_cis_window_file_reader(gene_cis_window_file)

                pheno_cov_numpy_dir = get_config()['associatr']['pheno_cov_numpy_dir']
                gene_pheno_cov = b.read_input(
                    f'{pheno_cov_numpy_dir}/{celltype}/{chromosome}/{gene}_pheno_cov.npy',
                )

                # Job 1: bcftools pre-filter to the SVs whose [POS, INFO/END]
                # intersects [cis_start, cis_end] by >=1 bp. Runs on the
                # bcftools image because the trtools image does not ship
                # bcftools/htslib. Output is a bgzip+tabix pair carried through
                # the batch as a resource group to the associaTR job below.
                filter_job = b.new_job(name=f'Filter SV VCF for {gene} [{celltype};{chromosome}]')
                if get_config()['associatr']['always_run']:
                    filter_job.always_run()
                filter_job.image(get_config()['images']['bcftools'])
                filter_job.storage(get_config()['associatr']['job_storage'])
                filter_job.cpu(get_config()['associatr']['job_cpu'])
                filter_job.declare_resource_group(
                    filtered_vcf={
                        'vcf.bgz': '{root}.vcf.bgz',
                        'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
                    },
                )
                filter_job.command(
                    f"bcftools view -r {cis_chrom}:1-{cis_end} "
                    f"-i 'INFO/END>={cis_start}' "
                    f"{variant_vcf.base} -Oz -o {filter_job.filtered_vcf['vcf.bgz']} && "
                    f"tabix -f -p vcf {filter_job.filtered_vcf['vcf.bgz']}",
                )

                # Job 2: associaTR on the per-gene filtered VCF. --region spans
                # the full chromosome so associaTR walks every surviving record.
                associatr_job = b.new_job(name=f'Run associatr on {gene} [{celltype};{chromosome}] SV')
                if get_config()['associatr']['always_run']:
                    associatr_job.always_run()
                associatr_job.image(get_config()['images']['trtools'])
                associatr_job.storage(get_config()['associatr']['job_storage'])
                associatr_job.cpu(get_config()['associatr']['job_cpu'])
                associatr_job.declare_resource_group(
                    association_results={'tsv': '{root}.tsv'},
                )
                associatr_job.command(
                    f"associaTR {associatr_job.association_results['tsv']} "
                    f"{filter_job.filtered_vcf['vcf.bgz']} "
                    f"{celltype}_{chromosome}_{gene} {gene_pheno_cov} "
                    f"--region={cis_chrom}:1-{chrom_len} --vcftype=eh",
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
    main()
