#!/usr/bin/env python3
"""
This script concatenates the results of running `coloc_runner.py` (output is per gene) into a single CSV file.

analysis-runner --dataset "bioheart" \
    --description "Parse coloc results" \
    --access-level "test" \
    --output-dir "str/associatr" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-bioheart-test/str/associatr/coloc \
    --celltypes=CD4_TCM \
    --phenos=ibd


"""
import click

from cpg_utils.hail_batch import get_batch


def coloc_results_combiner(coloc_dir, pheno, celltype):
    import pandas as pd

    from cpg_utils import to_path
    from cpg_utils.config import output_path

    files = list(to_path(f'{coloc_dir}/{pheno}/{celltype}').glob('*.tsv'))
    # List to store DataFrames
    dfs = []

    # Iterate over each file path
    for file_path in files:
        # Read file into a DataFrame
        df = pd.read_csv(file_path, sep='\t')
        # Append DataFrame to the list
        dfs.append(df)

    # Concatenate DataFrames row-wise
    result_df = pd.concat(dfs, ignore_index=True)
    # Write the result to a CSV file
    result_df.to_csv(output_path(f'coloc/{pheno}/{celltype}/gene_summary_result.csv', 'analysis'), index=False)


@click.command()
@click.option('--coloc-dir', help='Path to the directory containing coloc results')
@click.option('--celltypes', help='Cell type (can be multiple)')
@click.option('--phenos', help='Phenotype (can be multiple)')
def main(coloc_dir, celltypes, phenos):
    for celltype in celltypes.split(','):
        for pheno in phenos.split(','):
            b = get_batch()
            combiner_job = b.new_python_job(
                f'Coloc combiner for {celltype}:{pheno}',
            )
            combiner_job.call(coloc_results_combiner, coloc_dir, pheno, celltype)

    b.run(wait=False)


if __name__ == '__main__':
    main()
