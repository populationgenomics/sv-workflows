#!/usr/bin/env python3

import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path
import pandas as pd

"""
analysis-runner --dataset "bioheart" --description "Create a dataframe from all cell types" --access-level "test" \
    --output-dir "str/associatr/fine_mapping/susie_finemap" \
    --memory=highmem --cpu=16 --no-preemptible \
    pval_concatenator.py

"""

def main():

    cell_types = 'B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,CD4_TCM,gdT,dnT'

    cell_type_list = cell_types.split(',')
    dfs=[]
    for cell in cell_type_list:
        df = pd.read_csv(f'gs://gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/meta_results/min_p_0.05/{cell}/all_genes_pval_0.05.tsv', sep = '\t')
        dfs.append(df)
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/meta_results/min_p_0.05/all_genes_pval_0.05_noheader.tsv', sep = '\t', index = False, header=False)
    result_df.to_csv(f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/meta_results/min_p_0.05/all_genes_pval_0.05.tsv', sep = '\t', index = False)


if __name__ == '__main__':
    main()