#!/usr/bin/env python3

"""

This script will prepare the z and ld files for FINEMAP, based on outputs from associaTR (meta-analysis) and corr_matrix_maker.py.

The z file contains the effect size and standard error estimates for each variant associated with a gene.
The ld file contains the correlation matrix calcualtions for each variant associated with a gene.

analysis-runner --dataset "bioheart" \
    --description "Prepare files for FINEMAP" \
    --access-level "test" \
    --output-dir "str/associatr/fine_mapping" \
    finemap_files_prep.py \
    --ld-dir "gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/prep_files/v2-whole-copies-only/correlation_matrix" \
    --associatr-dir "gs://cpg-test-main-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results" \
    --celltypes "ASDC" \
    --chroms "chr22"



"""

import click
import numpy as np
import pandas as pd

import hailtop.batch as hb

from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch
from cpg_utils import to_path


def z_file_maker(gene_name, ld_file, associatr_dir, celltype, chrom):
    """
    Prepares the 'z file' for FINEMAP based on the output from associaTR (meta-analysis).

    """

    # Load the TSV file
    df = pd.read_csv(f'{associatr_dir}/{gene_name}_100000bp_meta_results.tsv', sep='\t')

    # Create the 'rsid' column
    df['rsid'] = df['chr'].astype(str) + ':' + df['pos'].astype(str) + '_' + df['motif']

    # Rename the 'chr' and 'pos' columns
    df.rename(columns={'chr': 'chromosome', 'pos': 'position', 'coeff_meta': 'beta', 'se_meta': 'se'}, inplace=True)

    # Create the 'allele1', 'allele2', 'maf' columns with NaN values
    df['allele1'] = 'nan'
    df['allele2'] = 'nan'
    df['maf'] = 'nan'

    # Select the required columns
    df = df[['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se']]

    # Load the index file
    df_index = pd.read_csv(ld_file, sep='\t')

    # Set the index of df to be 'rsid'
    df.set_index('rsid', inplace=True)

    # Reindex df according to the order of 'rsid' in df_index
    df = df.reindex(df_index['Unnamed: 0'])

    # Reset the index of df
    df.reset_index(inplace=True)
    df.rename(columns={'Unnamed: 0': 'rsid'}, inplace = True)

    output_file_path = output_path(f'finemap_prep/{celltype}/{chrom}/{gene_name}.z', 'analysis')
    # Save the DataFrame to a space-delimited file
    df.to_csv(output_file_path, sep=' ', index=False)

def ld_file_maker(gene_name, ld_file, celltype, chrom):
    """
    Prepares the 'ld file' for FINEMAP based on the output from corr_matrix_maker.py.
    """

    # Load the TSV file
    df = pd.read_csv(ld_file, sep='\t')
    df = df.drop(columns=['Unnamed: 0'])

    # Remove column names
    df.columns = ['' for _ in df.columns]

    # Remove index name
    df.index.names = [None]

    output_file_path = output_path(f'finemap_prep/{celltype}/{chrom}/{gene_name}.ld', 'analysis')
    # Save the DataFrame to a space-delimited file without index
    df.to_csv(output_file_path, sep=' ', index=False, header=False)

@click.option('--ld-dir', required=True, help='Path to the directory containing the LD files')
@click.option('--associatr-dir', required=True, help='Path to the directory containing the associatr files')
@click.option('--celltypes', required=True, help='Cell type')
@click.option('--chroms', required=True, help='Chromosome')
@click.option('--ld-job-cpu', help='CPU for LD job', default=0.25)
@click.option('--z-job-cpu', help='CPU for Z job', default=0.25)
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs to run', default=500)
@click.command()

def main(ld_dir, associatr_dir, celltypes, chroms, ld_job_cpu, z_job_cpu, max_parallel_jobs):

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    b = get_batch(name = f'FINEMAP prep for {celltypes}: {chroms}')
    for celltype in celltypes.split(','):
        for chrom in chroms.split(','):
            ld_files = list(to_path(f'{ld_dir}/{celltype}/{chrom}').glob('*.tsv'))
            for ld_file in ld_files:
                gene_name = str(ld_file).split('/')[-1].split('_')[0] #ENSG ID

                if to_path(
                    output_path(f'finemap_prep/{celltype}/{chrom}/{gene_name}.ld', 'analysis'),
                ).exists() and to_path(output_path(f'finemap_prep/{celltype}/{chrom}/{gene_name}.z', 'analysis')).exists():
                    continue
                ld_job = b.new_python_job(
                f'LD file maker for {celltype}:{chrom} {gene_name}',
            )
                ld_job.cpu(ld_job_cpu)
                ld_job.call(ld_file_maker, gene_name, ld_file, celltype, chrom)
                manage_concurrency_for_job(ld_job)
                z_job = b.new_python_job(f'Z file maker for {celltype}:{chrom} {gene_name}')
                z_job.cpu(z_job_cpu)
                z_job.call(z_file_maker, gene_name, ld_file, associatr_dir, celltype, chrom)
                manage_concurrency_for_job(z_job)

    b.run(wait=False)

if __name__ == '__main__':
    main()
