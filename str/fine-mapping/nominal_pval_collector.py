#!/usr/bin/env python3

import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path
import pandas as pd

"""
analysis-runner --dataset "bioheart" --description "Create a dataframe from all cell types" --access-level "test" \
    --output-dir "str/associatr/fine_mapping/susie_finemap" \
    nominal_pval_collector.py

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
            files_list = list(to_path(f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/meta_results/{cell}/chr{chrom}').glob('*.tsv'))
        except:
            continue
        for file in files_list:
            gene_name = str(file).split('/')[-1].split('_')[0]
            data = pd.read_csv(file, sep = '\t')
            data['celltype'] = cell
            data['gene'] = gene_name
            data = data[data['pval_meta'] < 0.05] # Filter out non-significant results (0.05)
            data['motif_len'] = data['motif'].str.len()
            data['end'] = (data['pos'].astype(float) + data['ref_len'].astype(float) * data['motif_len'].astype(float)).round().astype(int)
            data = data[['chr', 'pos','end' ,'motif', 'celltype','coeff_meta', 'se_meta', 'pval_meta']]
            dfs.append(data)
        print (f'Processed {cell} and chr{chrom}')
    result_df = pd.concat(dfs, ignore_index=True)
    result_df.to_csv(f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/meta_results/min_p_0.05/{cell}/all_genes_pval_0.05.tsv', sep = '\t', index = False)

def main():

    cell_types = 'B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,CD4_TCM,gdT,dnT'
    b = get_batch()
    # Split the string into a list of characters
    cell_type_list = cell_types.split(',')
    for cell in cell_type_list:
        combiner_job = b.new_python_job(
            f'Concatenate all genes for {cell}',
        )
        combiner_job.cpu(2)
        combiner_job.call(run_concatenator, cell)
    b.run(wait=False)



if __name__ == '__main__':
    main()