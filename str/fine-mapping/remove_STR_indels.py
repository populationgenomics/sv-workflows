#!/usr/bin/env python3

"""
This script filters raw associaTR outputs that contain both eSTR and eSNP results (i.e. we assume that `dataframe_concatenator.py` has been run).

Particularly, we remove indels that actually represent STRs, AND remove duplicate eSTRs (retain only one eSTR per duplicate set).
This is necessary to improve the accuracy of fine-mapping.

Indels are considered STRs (and subsequently removed) if:

- they overlap with an STR region, specified by an eSTR in the associaTR output; and
- the (reverse complement) sequence of the indel is a whole copy of at least one cyclical representation of the STR motif.
e.g. a 'TAT' insertion overlapping with a 'ATT' STR would be considered an STR as 'TAT' is a cyclical representation of 'ATT'.
 a 'GC' insertion overlapping with a 'CAG' STR would not be considered an STR, as 'GC' is a partial copy of 'GCA', a cyclical representation of 'CAG'.

Note that impure indels are conservatively not considered STRs. For example, a 'GCCGCA' insertion overlapping a 'GCC' STR would not be considered an STR.

This script additionally removes duplicate eSTRs (defined by sharing the same coordinates and motif), retaining only one eSTR per duplicate set (chosen based on having the lowest p-value).

Usage:

analysis-runner --dataset bioheart --output-dir str/associatr/snps_and_strs/rm_str_indels_dup_strs/tob_n1055_and_bioheart_n990/meta_results \
--description "Remove STR indels and duplicate eSTRs" \
--access-level "test" \
remove_STR_indels.py \
--associatr-dir=gs://cpg-bioheart-test/str/associatr/snps_and_strs/tob_n1055_and_bioheart_n990/meta_results/meta_results \
--celltypes=gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive \
--chromosomes=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22




"""
import ast

import click
import numpy as np

import hailtop.batch as hb
from hailtop.batch import ResourceGroup

from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch


def check_str(motif):
    """
    Checks if a row in a df represents an STR (i.e. the motif does not contain a '-')
    """
    return '-' not in motif


def reverse_complement(sequence):
    """
    Return the reverse complement of a DNA sequence.

    """
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
    return reverse_complement_sequence


def cyclical_shifts(s):
    """
    Generate all cyclical shifts of a string.
    """
    return [s[i:] + s[:i] for i in range(len(s))]


def filter_str_indels_and_duplicates(associatr_dir, celltype, chrom):
    import pandas as pd

    from cpg_utils import to_path

    # Load associaTR output
    gene_files = list(to_path(f'{associatr_dir}/{celltype}/{chrom}').glob('*.tsv'))
    for gene_file in gene_files:
        data = pd.read_csv(str(gene_file), sep='\t')
        gene_file_name = str(gene_file).split('/')[-1]
        print(f'Processing {gene_file_name}...')

        # Find duplicate eSTRs and remove all but the one with the lowest p-value
        duplicates = data[data.duplicated(subset=['chr', 'pos', 'motif'], keep=False)]
        # Group by 'chr', 'pos', 'motif' and keep the one with the lowest 'p-val'
        lowest_pval_duplicates = duplicates.loc[duplicates.groupby(['chr', 'pos', 'motif'])['pval_meta_fixed'].idxmin()]
        # Find non-duplicate rows
        non_duplicates = data.drop(duplicates.index)
        # Concatenate the non-duplicates with the lowest p-value duplicates
        result_df = pd.concat([non_duplicates, lowest_pval_duplicates]).sort_index()

        # get df containing only STRs
        result_df_str = result_df[result_df['motif'].apply(check_str)]
        result_df_str.loc[:, 'end'] = (
            result_df_str.loc[:, 'pos'] + result_df_str.loc[:, 'ref_len'] * result_df_str.loc[:, 'period']
        )

        indices_to_delete = []  # store indices of indels representing STRs (to drop later)
        for i, row1 in result_df.iterrows():  # iterate through every variant record in df
            if '-' not in row1['motif']:  # skip if STR
                continue
            parts = row1['motif'].split('-')
            if len(parts[0]) == 1 and len(parts[1]) == 1:  # skip if SNP (single base change)
                continue

            pos_indel = row1['pos']
            # Split the motif to find the indel
            motif_parts = row1['motif'].split('-')

            indel = next((part for part in motif_parts if len(part) != 1), None)[1:]

            found = False
            for j, row2 in result_df_str.iterrows():
                if (
                    row2['pos'] <= pos_indel <= row2['end']
                ):  # Check if pos_indel is between 'pos' and 'end' of any row in STR dataframe
                    str_motif = row2['motif']

                    if (indel == str_motif) or (
                        reverse_complement(indel) == str_motif
                    ):  # straightforward case where (reverse complement of) indel matches STR motif exactly
                        # print(f"Indel: {row1['motif']} {indel} matches motif {str_motif}")
                        indices_to_delete.append(i)
                        found = True
                        break

                    if len(indel) % len(str_motif) == 0:  # check if len(indel) is a multiple of len(str_motif)
                        elongated_motif = str_motif * (len(indel) // len(str_motif))
                        if (indel in cyclical_shifts(elongated_motif)) or (
                            reverse_complement(indel) in cyclical_shifts(elongated_motif)
                        ):
                            indices_to_delete.append(i)
                            found = True
                            break

                    found = True
                    break
            # if not found:
            # print("Indel not in STR interval")
        # Drop the rows where indels are representing STRs
        result_df = result_df.drop(indices_to_delete)
        result_df.to_csv(output_path(f'{celltype}/{chrom}/{gene_file_name}', 'analysis'), sep='\t', index=False)


@click.option('--associatr-dir', required=True, type=str, help='Directory containing associaTR outputs.')
@click.option('--celltypes', required=True, type=str, help='Comma-separated list of cell types to process.')
@click.option('--chromosomes', required=True, type=str, help='Comma-separated list of chromosomes to process.')
@click.option(
    '--max-parallel-jobs',
    required=True,
    type=int,
    help='Maximum number of jobs to run in parallel.',
    default=50,
)
@click.option('--job-cpu', required=False, type=float, help='Number of CPUs to use per job.', default=0.25)
@click.option('--job-storage', required=False, type=str, help='Storage to use per job.', default='0G')
@click.command()
def main(
    associatr_dir: str,
    celltypes: str,
    chromosomes: str,
    max_parallel_jobs: int,
    job_cpu: int,
    job_storage: str,
):
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    b = get_batch(name='Remove STR indels and duplicate eSTRs')
    for celltype in celltypes.split(','):
        for chrom in chromosomes.split(','):
            filter_job = b.new_python_job(
                f'Remove STR-indels and duplicate eSTRs {celltype}:{chrom}',
            )
            filter_job.cpu(job_cpu)
            filter_job.storage(job_storage)

            filter_job.call(filter_str_indels_and_duplicates, associatr_dir, celltype, chrom)
            manage_concurrency_for_job(filter_job)

    b.run(wait=False)


if __name__ == '__main__':
    main()
