#!/usr/bin/env python3
"""
This script concatenates the results of running `coloc_runner.py` (output is per gene) into a single CSV file.

analysis-runner --dataset "bioheart" \
    --description "Parse coloc results" \
    --access-level "test" \
    --memory='16G' \
    --storage='10G' \
    --output-dir "str/associatr" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-bioheart-test/str/associatr/freeze_1/gwas_ld/bioheart-only-snps/gwas_ld \
    --celltypes=gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive \
    --phenos=ibd


"""
import click

from cpg_utils.hail_batch import get_batch


def coloc_results_combiner(coloc_dir, pheno, celltype):
    import pandas as pd

    from cpg_utils import to_path
    from cpg_utils.config import output_path

    files = list(to_path(f'{coloc_dir}/{pheno}/{celltype}').glob('*.csv'))
    # List to store DataFrames
    dfs = []

    # Iterate over each file path
    for file_path in files:
        # Read file into a DataFrame
        df = pd.read_csv(file_path)
        # Append DataFrame to the list
        dfs.append(df)

    # Concatenate DataFrames row-wise
    result_df = pd.concat(dfs, ignore_index=True)
    # Write the result to a CSV file
    result_df.to_csv(f'{coloc_dir}/{pheno}/{celltype}/gene_summary_result.csv', index=False)


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
