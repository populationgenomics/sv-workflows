#!/usr/bin/env python3
"""
This script concatenates the results of running `coloc_runner.py` (output is per gene) into a single CSV file.

analysis-runner --dataset "bioheart" \
    --description "Parse coloc results" \
    --access-level "test" \
    --output-dir "tenk10k/str/associatr/final_freeze/coloc-tr-snp/sig_str_and_gwas_hit" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/coloc-tr-snp/sig_str_and_gwas_hit \
    --celltypes=B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,gdT,dnT,CD4_TCM \
    --phenos=Trujillo_methylation_eQTLs


    cpg-bioheart-test-analysis/str/associatr/coloc-snp-only/sig_str_and_gwas_hit/lymphocytic_leukemia_GCST90011814

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
        df['probe'] = str(file_path).split('/')[-1].split('_')[1]
        # Append DataFrame to the list
        dfs.append(df)

    # Concatenate DataFrames row-wise
    result_df = pd.concat(dfs, ignore_index=True)
    # Write the result to a CSV file
    result_df.to_csv(output_path(f'{pheno}/{celltype}/gene_summary_result.csv', 'analysis'), index=False)


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
