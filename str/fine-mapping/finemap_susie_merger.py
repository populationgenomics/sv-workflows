#!/usr/bin/env python3
"""
This script merges the output of FINEMAP and SUSIE into a single file per gene-celltype combination.

analysis-runner --dataset "bioheart" --access-level 'test' --description "Merge FINEMAP and SUSIE results" \
    --output-dir "str/associatr/fine_mapping" \
    finemap_susie_merger.py \
    --finemap-dir "gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/finemap/ofiles" \
    --susie-dir "gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/v2-whole-copies-only-v2/susie" \
    --celltypes "B_memory" \
    --chromosomes "chr4"


"""


import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path


def run_concatenator(finemap_dir, susie_dir, celltype, chromosome, gene):
    """
    Concatenate two dataframes together.
    """
    import pandas as pd

    from cpg_utils.hail_batch import output_path

    # read input files
    finemap_df = pd.read_csv(f'{finemap_dir}/{celltype}/{chromosome}/{gene}.snp', sep='\t')
    finemap_df = finemap_df[['rsid', 'beta', 'se', 'prob', 'log10bf']]
    finemap_df = finemap_df.rename(columns={'prob': 'finemap_prob', 'log10bf': 'finemap_log10bf'})
    susie_df = pd.read_csv(f'{susie_dir}/{celltype}/{chromosome}/{gene}_100kb.tsv', sep='\t')

    # Create new columns with rounded values
    susie_df['beta_rounded'] = susie_df['coeff_meta'].round(3)
    susie_df['se_rounded'] = susie_df['se_meta'].round(3)

    finemap_df['beta_rounded'] = finemap_df['beta'].round(3)
    finemap_df['se_rounded'] = finemap_df['se'].round(3)

    # Merge on the new columns
    merged_df = susie_df.merge(
        finemap_df,
        left_on=['varid', 'beta_rounded', 'se_rounded'],
        right_on=['rsid', 'beta_rounded', 'se_rounded'],
    )

    # drop the new columns as they are no longer needed
    merged_df = merged_df.drop(columns=['beta_rounded', 'se_rounded'])

    # write results as a tsv file to gcp
    merged_df.to_csv(
        output_path(f'susie_finemap/{celltype}/{chromosome}/{gene}.tsv', 'analysis'), sep='\t', index=False,
    )


@click.option('--finemap-dir', help='Input directory for the FINEMAP .snp files')
@click.option('--susie-dir', help='Input directory for the susie output .tsv files')
@click.option('--celltypes', help='comma-separated list of cell types')
@click.option('--chromosomes', help='comma-separated list of chromosomes')
@click.option('--max-parallel-jobs', help='Maximum number of jobs to run in parallel', default=500)
@click.option('--job-cpu', help='Number of CPUs to use for each job', default=0.25)
@click.option('--always-run', help='Job set to always run', is_flag=True)
@click.command()
def main(finemap_dir, susie_dir, celltypes, chromosomes, max_parallel_jobs, always_run, job_cpu):
    """
    Runner script to FINEMAP and SUSIE DFs together
    """
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    b = get_batch(name='Merge FINEMAP and SUSIE results')
    for celltype in celltypes.split(','):
        for chromosome in chromosomes.split(','):
            # gets the list of available gene files in input dir 1
            gene_files = list(map(str, to_path(f'{susie_dir}/{celltype}/{chromosome}').glob('*.tsv')))
            for gene_file in gene_files:
                # see if the gene file exists in input dir 2. If it does, concatenate the two files together.
                gene = gene_file.split('/')[-1].split('_')[0]
                finemap_file_name = gene + '.snp'
                if to_path(f'{finemap_dir}/{celltype}/{chromosome}/{finemap_file_name}').exists():
                    # see if output file exists. If it does, skip this job
                    if to_path(output_path(f'susie_finemap/{celltype}/{chromosome}/{gene}.tsv', 'analysis')).exists():
                        continue
                    j = b.new_python_job(name=f'Merge SUSIE and FINEMAP results for {celltype}_{chromosome}_{gene}')
                    j.cpu(job_cpu)
                    if always_run:
                        j.always_run()
                    j.call(run_concatenator, finemap_dir, susie_dir, celltype, chromosome, gene)
                    manage_concurrency_for_job(j)
    b.run(wait=False)


if __name__ == '__main__':
    main()
