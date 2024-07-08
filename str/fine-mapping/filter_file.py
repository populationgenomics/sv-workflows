#!/usr/bin/env python3
"""
This script concatenates the results of running `coloc_runner.py` (output is per gene) into a single CSV file.

analysis-runner --dataset "bioheart" \
    --description "Parse coloc results" \
    --access-level "test" \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr" \
    filter_file.py \
    --celltypes "B_intermediate"
"""
import click
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def fold_calc(directory, celltype):
    import pandas as pd

    from cpg_utils import to_path
    from cpg_utils.config import output_path

    gene_files = list(to_path(f'{directory}/{celltype}').rglob('*.tsv'))

    # List to store DataFrames
    dfs = []

    # Iterate over each file path
    for file_path in gene_files:
        gene_name = file_path.split('/')[-1].split('_')[0]
        # Read file into a DataFrame
        df = pd.read_csv(file_path, sep='\t')
        str_df = df[~df['motif'].str.contains('-')]
        snp_df = df[df['motif'].str.contains('-')]
        str_df_min = str_df['pval_meta'].min()
        snp_df_min = snp_df['pval_meta'].min()
        pval_fold = str_df_min / snp_df_min
        result_df = pd.DataFrame(
            {
                'gene': [gene_name],
                'str_lead_pval': [str_df_min],
                'snp_lead_pval': [snp_df_min],
                'pval_fold': [pval_fold],
                'celltype': celltype,
            },
        )
        # Append DataFrame to the list
        dfs.append(result_df)

    # Concatenate DataFrames row-wise
    master_result_df = pd.concat(dfs, ignore_index=True)
    # Write the result to a CSV file
    master_result_df.to_csv(
        f'gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results/pval_fold_enrichment/{celltype}_results.csv',
        index=False,
    )

@click.option('--celltypes', help='Cell type (can be multiple)')
@click.command()
def main(celltypes):
    b = get_batch()
    for celltype in celltypes.split(','):
        directory = 'gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results'
        combiner_job = b.new_python_job(
            f'Parser job for {celltype}',
        )
        combiner_job.cpu(4)
        combiner_job.call(fold_calc, directory, celltype)

    b.run(wait=False)


if __name__ == '__main__':
    main()
