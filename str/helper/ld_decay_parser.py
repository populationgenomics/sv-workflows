#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script writes out genes for every cell type if the gene has a lead signal that is not a SNP.
We test all genes tested, not just eGenes that pass an FDR.
Used to generate one stat in the paper.

 analysis-runner  --dataset "bioheart" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr/estrs" \
python3 ld_decay_parser.py

"""
import json

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

def cell_parser(cell, chrom_num):
    dfs = []
    files =list(to_path(f'gs://cpg-bioheart-test-analysis/str/associatr/estrs/ld_decay/test/mut_ex/{cell}/chr{chrom_num}').rglob(f'{cell}_chr{chrom_num}_ENSG*'))
    for file in files:
        try:
            df = pd.read_csv(file, sep = '\t')

        except:
            continue

        dfs.append(df)

    result_df = pd.concat(dfs, ignore_index = True)
    result_df.to_csv(f'gs://cpg-bioheart-test-analysis/str/associatr/estrs/ld_decay/test/mut_ex/{cell}/chr{chrom_num}/all_genes_ld_decay_metrics.csv', index = False)



def main():
    cell_types = 'B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,gdT,dnT,CD4_TCM'
    b = get_batch(name='LD decay parser')
    # Split the string into a list of characters
    cell_type_list = cell_types.split(',')
    for cell in cell_type_list:
        for chrom_num in [22]:
            j = b.new_python_job(
                name=f'Get LD decay metrics for {cell} and chr{chrom_num}',
            )
            j.call(
                cell_parser,
                cell,
                chrom_num,
            )

    #result_df = pd.concat(dfs, ignore_index = True)
    #result_df.to_csv('gs://cpg-bioheart-test-analysis/str/associatr/estrs/ld_decay/test/mut_ex/chr22_all_cell_types_ld_decay_metrics.csv', index = False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments
