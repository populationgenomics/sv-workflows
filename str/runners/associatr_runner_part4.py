#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script is step 3 of 4 for running associaTR.
It aims to:
- filter VCF to intersect with cis window of a gene
- run associatr on the Gene + cell type with STR genotypes

 analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files for associatr" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
     associatr_runner_part2.py  --celltypes=Plasma --chromosomes=chr22

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



# inputs:
@click.option('--celltypes')
@click.option('--chromosomes', help=' eg chr22')
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=50,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
@click.command()
def main(
    celltypes, chromosomes, max_parallel_jobs, cis_window_size
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
            variant_vcf = b.read_input("gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/dumpSTR/filtered_mergeSTR_results.vcf")
            with to_path(output_path(f'input_files/scRNA_gene_lists/{celltype}/{chromosome}_{celltype}_filtered_genes.json')).open('r') as file:
                pseudobulk_gene_names = json.load(file)
            for gene in pseudobulk_gene_names:
                # intersect VCF with cis window of the gene
                cis_window_file = b.read_input(f'gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/cis_window_files/{celltype}/{chromosome}/{gene}_100000bp.bed')
                #need to extract the gene start and end from the cis window file for input into 'region'

                #run associaTR job on the gene
                #associatr_job = b.new_job(name=f'Run associatr on {gene} [{celltype};{chromosome}]')
                #associatr_job.image(get_config()['images']['trtools'])
                #associatr_job.storage('8G')
                #associatr_job.depends_on(bedtools_job)
                #associatr_job.command(
                #    f' associaTR association_results.tsv {bedtools_job.ofile} {celltype}_{chromosome}_{gene}  '
                #)

                b.write_output(gene_cis_job.ofile, output_path(f'input_files/cis_window_files/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.bed'))
                manage_concurrency_for_job(gene_cis_job)
    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter






#dumpSTR --vcf mergeSTR_1057_samples_eh.vcf.gz --out filtered_mergeSTR_results --vcftype eh --min-locus-callrate 0.9 --min-locus-het 0.1 --min-locus-hwep 0.0001 --filter-regions segDupRegions_hg38_sorted.bed.gz
