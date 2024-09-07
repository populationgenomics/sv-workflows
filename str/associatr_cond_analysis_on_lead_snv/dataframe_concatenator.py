#!/usr/bin/env python3

"""
This script concatenates (row-wise) two dataframes together.
Assumed that the header for the two dataframes are the same.
Outputs results as a TSV file.

This helper script will be used to concatenate the results of the meta-analysis for eSNPs and eSTRs together, with each gene having its own file.
Only common genes between the two datasets will be concatenated.

analysis-runner --dataset "bioheart" --description "concatenate meta-analysis results" --access-level "test" \
    --output-dir "str/associatr/cond_analysis_lead_snv/snps_and_strs/tob_n1055_and_bioheart_n990/meta_results" \
    dataframe_concatenator.py \
    --input-dir-1=gs://cpg-bioheart-test-analysis/str/associatr/cond_analysis_lead_snv/common_variants_snps/tob_n1055_and_bioheart_n990/meta_results \
    --input-dir-2=gs://cpg-bioheart-test-analysis/str/associatr/cond_analysis_lead_snv/tob_n1055_and_bioheart_n990/meta_results \
    --celltypes=gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive

"""

import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path


def run_concatenator(input_dir_1, input_dir_2, celltype, chromosome, gene_file):
    """
    Concatenate two dataframes together.
    """
    import pandas as pd

    from cpg_utils.hail_batch import output_path

    # read input files
    df1 = pd.read_csv(f'{input_dir_1}/{celltype}/{chromosome}/{gene_file}', sep='\t')
    df2 = pd.read_csv(f'{input_dir_2}/{celltype}/{chromosome}/{gene_file}', sep='\t')
    # concatenate
    df = pd.concat([df1, df2], ignore_index=True)
    # write results as a tsv file to gcp
    df.to_csv(output_path(f'{celltype}/{chromosome}/{gene_file}', 'analysis'), sep='\t', index=False)


@click.option('--input-dir-1', help='Input directory for the first collection of dataframes')
@click.option('--input-dir-2', help='Input directory for the second collection of dataframes')
@click.option('--celltypes', help='comma-separated list of cell types')
@click.option('--max-parallel-jobs', help='Maximum number of jobs to run in parallel', default=500)
@click.option('--always-run', help='Job set to always run', is_flag=True)
@click.command()
def main(input_dir_1, input_dir_2, celltypes, max_parallel_jobs, always_run):
    """
    Runner script to concatenate two dataframes together.
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

    b = get_batch()
    for celltype in celltypes.split(','):
        # gets the list of available gene files in input dir 1
        gene_files = list(map(str, to_path(f'{input_dir_1}/{celltype}').rglob('*.tsv')))
        for gene_file in gene_files:
            # see if the gene file exists in input dir 2. If it does, concatenate the two files together.
            file_name = gene_file.split('/')[-1]
            chromosome = gene_file.split('/')[-2]
            if to_path(f'{input_dir_2}/{celltype}/{chromosome}/{file_name}').exists():
                # see if output file exists. If it does, skip this job
                if to_path(output_path(f'{celltype}/{chromosome}/{file_name}', 'analysis')).exists():
                    continue
                j = b.new_python_job(name=f'concatenate_{celltype}_{chromosome}_{file_name}')
                j.cpu(0.25)
                if always_run:
                    j.always_run()
                j.call(run_concatenator, input_dir_1, input_dir_2, celltype, chromosome, file_name)
                manage_concurrency_for_job(j)
    b.run(wait=False)


if __name__ == '__main__':
    main()
