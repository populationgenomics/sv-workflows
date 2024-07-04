#!/usr/bin/env python3
"""
This script concatenates the results of running `coloc_runner.py` (output is per gene) into a single CSV file.

analysis-runner --dataset "bioheart" \
    --description "Parse coloc results" \
    --access-level "test" \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr" \
    coloc_results_parser.py \
    --coloc-dir=gs://cpg-bioheart-test-analysis/str/associatr/coloc/sig_str_and_gwas_hit \
    --celltypes=gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive \
    --phenos=alanine-aminotransferase,albumin,alkaline_phosphatase,apolipoprotein_a,apolipoprotein_b,aspartate_aminotransferase,c_reactive_protein,calcium,cholesterol,creatinine,cystatin_c,eosinophil_count,eosinophil_percent,gamma_glutamyltransferase,glucose,haematocrit,haemoglobin_concentration,lymphocyte_count,lymphocyte_percent,mean_corpuscular_haemoglobin,mean_corpuscular_haemoglobin_concentration,mean_corpuscular_volume,mean_platelet_volume,mean_sphered_cell_volume,neutrophil_count,neutrophil_percent,platelet_count,platelet_crit,platelet_distribution_width,red_blood_cell_count,red_blood_cell_distribution_width,white_blood_cell_count


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
    result_df.to_csv(output_path(f'coloc-snp/sig_str_and_gwas_hit/{pheno}/{celltype}/gene_summary_result.csv', 'analysis'), index=False)


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
