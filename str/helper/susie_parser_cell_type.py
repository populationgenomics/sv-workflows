#!/usr/bin/env python3

import pandas as pd
from cpg_utils import to_path, get_batch

cell_types = 'B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,CD4_TCM,gdT,dnT'

# Split the string into a list of cell types
cell_type_list = cell_types.split(',')

b = get_batch(name='SuSiE Parser')

for cell in cell_type_list:

    # Define the job
    susie_job = b.new_python_job(f'Parse SuSiE results for {cell}')
    susie_job.storage('50G')
    susie_job.memory('32G')

    # Define the function to be called by the job
    def parse_susie(cell):
        files_list = list(to_path(f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/fine_mapping/susie/{cell}').rglob('*.tsv'))
        dfs = []
        for file in files_list:
            gene_name = str(file).split('/')[-1].split('_')[0]
            data = pd.read_csv(file, sep='\t')
            data['celltype'] = cell
            data['gene'] = gene_name
            dfs.append(data)
        result_df = pd.concat(dfs, ignore_index=True)
        result_df.to_csv(f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/fine_mapping/susie/{cell}/all_genes.csv', index=False)

    # Call the function within the job
    susie_job.call(parse_susie, cell)

    # Submit the batch
b.run(wait=False)