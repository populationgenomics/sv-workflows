#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script runs associaTR, given a chromosome or cell type.
Ensure prior scripts have been run to generate dependent files, particularly:
- get_cis_numpy_files.py
- pseudobulk.py
- get_covariates.py
- qc_filters_associatr.py (depends on qc_annotator.py)

 analysis-runner --dataset "bioheart" \
    --description "run associatr" \
    --access-level "test" \
    --output-dir "str/associatr" \
     associatr_runner_part5.py  --celltypes=ASDC --chromosomes=chr1 \
    --vcf-file-path=gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1/bctools/hail_filtered.vcf.gz
    --gene-list-dir=gs://cpg-bioheart-test/str/associatr/input_files/scRNA_gene_lists/0.01_non_zero_threshold
    --cis-window-dir=gs://cpg-bioheart-test/str/associatr/input_files/cis_window_files/v1
    --pheno-cov-numpy-dir=gs://cpg-bioheart-test/str/associatr/input_files/pheno_cov_numpy/v1


"""
import json
import click
import pandas as pd
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path, init_batch


def gene_cis_window_file_reader(file_path):
    cis_window = pd.read_csv(file_path, sep='\t', header=None)
    cis_window.columns = ['chrom', 'start', 'end']
    return f'{cis_window["chrom"][0]}:{cis_window["start"][0]}-{cis_window["end"][0]}'


# inputs:
@click.option('--celltypes')
@click.option('--chromosomes', help=' eg chr22')
@click.option(
    '--vcf-file-path', help='gs://... to the output of qc_filters_associatr.py'
)
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=50,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
@click.option('--cis-window-size', type=int, default=100000)
@click.option('--gene-list-dir', help='directory to gene list json files')
@click.option('cis-window-dir', help='directory to cis window files')
@click.option('--pheno-cov-numpy-dir', help='directory to numpy files')
@click.option('--version', help='version of the output', default='v1')
@click.command()
def main(
    celltypes,
    chromosomes,
    max_parallel_jobs,
    cis_window_size,
    vcf_file_path,
    gene_list_dir,
    cis_window_dir,
    pheno_cov_numpy_dir,
    version,
):
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
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    for celltype in celltypes.split(','):
        for chromosome in chromosomes.split(','):
            variant_vcf = b.read_input_group(
                **dict(
                    base=vcf_file_path,
                    tbi=vcf_file_path + '.tbi',
                )
            )
            with to_path(
                output_path(
                    f'{gene_list_dir}/{celltype}/{chromosome}_{celltype}_gene_list.json'
                )
            ).open('r') as file:
                pseudobulk_gene_names = json.load(file)
            for gene in pseudobulk_gene_names:
                gene_cis_window_file = f'{cis_window_dir}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.bed'
                # need to extract the gene start and end from the cis window file for input into 'region'
                cis_window_region = gene_cis_window_file_reader(gene_cis_window_file)
                gene_pheno_cov = b.read_input(
                    f'{pheno_cov_numpy_dir}/{celltype}/{chromosome}_{celltype}_{gene}.npy'
                )

                # run associaTR job on the gene
                associatr_job = b.new_job(
                    name=f'Run associatr on {gene} [{celltype};{chromosome}]'
                )
                associatr_job.image(get_config()['images']['trtools'])
                associatr_job.storage('8G')
                associatr_job.cpu(4)
                associatr_job.declare_resource_group(
                    association_results={'tsv': '{root}.tsv'}
                )
                associatr_job.command(
                    f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.base} {celltype}_{chromosome}_{gene} {gene_pheno_cov} --region={cis_window_region} --vcftype=eh"
                )

                b.write_output(
                    associatr_job.association_results,
                    output_path(
                        f'results/{version}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp'
                    ),
                    'analysis',
                )
                manage_concurrency_for_job(associatr_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
