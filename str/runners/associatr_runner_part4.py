#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script is step 4 of 4 for running associaTR.
It aims to:
- run associatr on the Gene + cell type with STR genotypes

 analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files for associatr" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
     associatr_runner_part4.py  --celltypes=Plasma --chromosomes=chr22 \
    --vcf-file-path=gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/hail/hail_filtered.vcf.gz

"""
import click
import pandas as pd
import hail as hl
import hailtop.batch as hb
import json

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch


from cpg_utils.hail_batch import output_path, init_batch

config = get_config()

def gene_cis_window_file_reader(file_path):
    cis_window = pd.read_csv(file_path, sep='\t', header=None)
    cis_window.columns = ['chrom', 'start', 'end']
    return f'{cis_window["chrom"][0]}:{cis_window["start"][0]}-{cis_window["end"][0]}'

# inputs:
@click.option('--celltypes')
@click.option('--chromosomes', help=' eg chr22')
@click.option('--vcf-file-path', help='gs://... to the output of dumpSTR')
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=50,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
@click.option('--cis-window-size', type=int,default=100000)
@click.command()
def main(
    celltypes, chromosomes, max_parallel_jobs, cis_window_size,vcf_file_path
):
    """
    Run associaTR processing pipeline
    """
    config = get_config()
    b = get_batch()
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
            variant_vcf =b.read_input_group(
            **dict(
            base = vcf_file_path,
            tbi = vcf_file_path + '.tbi',
            )
            )
            with to_path(output_path(f'input_files/scRNA_gene_lists/{celltype}/{chromosome}_{celltype}_filtered_genes.json')).open('r') as file:
                pseudobulk_gene_names = json.load(file)
            for gene in ['MYH9', 'RTCB', 'TSPO', 'ADSL','PES1']:
                gene_cis_window_file = f'gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/cis_window_files/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.bed'
                cis_window_region = gene_cis_window_file_reader(gene_cis_window_file)
                gene_pheno_cov = b.read_input(f'gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/gene_pheno_cov_numpy/{celltype}/{chromosome}_{celltype}_{gene}.npy')
                #need to extract the gene start and end from the cis window file for input into 'region'

                #run associaTR job on the gene
                associatr_job = b.new_job(name=f'Run associatr on {gene} [{celltype};{chromosome}]')
                associatr_job.image(get_config()['images']['trtools'])
                associatr_job.storage('8G')
                associatr_job.cpu(4)
                associatr_job.declare_resource_group(
                    association_results={'tsv': '{root}.tsv'}
                )
                associatr_job.command(
                    f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.base} {celltype}_{chromosome}_{gene} {gene_pheno_cov} --region={cis_window_region} --vcftype=eh"
                )

                b.write_output(associatr_job.association_results, output_path(f'output-files/{celltype}/{chromosome}/{gene}_{cis_window_size}bp'))
                manage_concurrency_for_job(associatr_job)
    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter






#dumpSTR --vcf mergeSTR_1057_samples_eh.vcf.gz --out filtered_mergeSTR_results --vcftype eh --min-locus-callrate 0.9 --min-locus-het 0.1 --min-locus-hwep 0.0001 --filter-regions segDupRegions_hg38_sorted.bed.gz
