#!/usr/bin/env python3

import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path
import pandas as pd

"""
analysis-runner --dataset "bioheart" --description "Create a dataframe from all cell types" --access-level "test" \
    --output-dir "str/associatr/fine_mapping/susie_finemap" \
    --memory=16G --storage=30G \
    df_maker.py

"""

def run_concatenator(cell):
    """
    Concatenate two dataframes together.
    """
    import pandas as pd

    # read input files
    dfs =[]
    for chrom in range(1,23):
        try:
            files_list = list(to_path(f'gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/{cell}/chr{chrom}').glob('*.tsv'))
        except:
            continue
        for file in files_list:
            gene_name = str(file).split('/')[-1].split('_')[0]
            data = pd.read_csv(file, sep = '\t')
            data['celltype'] = cell
            data['gene'] = gene_name
            dfs.append(data)
        print (f'Processed {cell} and chr{chrom}')
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/{cell}/all_genes.tsv', sep = '\t', index = False)

def main():

    cell_types = 'B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,CD4_TCM,gdT,dnT'

    # Split the string into a list of characters
    cell_type_list = cell_types.split(',')
    dfs = []
    for cell in cell_type_list:
        data = pd.read_csv(f'gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/{cell}/all_genes.tsv', sep = '\t')
        dfs.append(data)
    result_df = pd.concat(dfs, ignore_index=True)
    result_df = result_df[result_df['pval_meta'] < 5e-8]
    result_df.to_csv(f'gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/all_cell_types_all_genes_sig_only.tsv', sep = '\t', index = False)



if __name__ == '__main__':
    main()