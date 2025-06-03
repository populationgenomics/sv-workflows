#!/usr/bin/env python3
"""
This script concatenates the results from running SUSIE-R on a cell type basis.
The output is a single CSV file containing per cell type results (optional min PIP filtering)

analysis-runner --dataset bioheart --access-level test --output-dir tenk10k/str/associatr/final_freeze/fine_mapping/susie_mcv/output_files --description "Concatenate SuSIE results" python3 results_concatenator.py \
    --celltypes "CD14_Mono,ASDC,cDC1,gdT,NK,B_naive,B_memory,CD4_CTL,CD4_TCM,CD8_TCM,cDC2,HSPC,NK_CD56bright,Plasmablast,B_intermediate,CD4_Naive,CD4_TEM,CD8_TEM,dnT,MAIT,NK_Proliferating,Treg,CD16_Mono,CD4_Proliferating,CD8_Naive,pDC,CD8_Proliferating,ILC"


"""

import click
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch


def process_file(file, pip_prune_threshold):
    import pandas as pd

    df = pd.read_csv(file, sep='\t')
    cell_type = str(file).split('/')[-2]
    df['cell_type'] = cell_type  # add a column for the cell type
    df['gene_name']= file.split('/')[-1].split('_')[0]  # extract gene name from file name
    df = df[df['pip_prune'] >= pip_prune_threshold]  # actually filter
    df.drop(columns=['chr'], inplace=True)  # drop the chr column
    return df

def concatenator(input_dir, cell_type, pip_prune_threshold):
    import pandas as pd
    from concurrent.futures import ThreadPoolExecutor
    from cpg_utils.hail_batch import output_path

    files = list(to_path(f'{input_dir}/{cell_type}').glob('*.tsv'))
    master_df = pd.DataFrame()
    with ThreadPoolExecutor() as executor:
        # Pass pip_prune_threshold to process_file
        results = list(executor.map(lambda f: process_file(f, pip_prune_threshold), files))
    master_df = pd.concat(results, ignore_index=True)
    output_gcs = output_path(f'concat_results/{cell_type}_concatenated_results.tsv', 'analysis')
    master_df.to_csv(output_gcs, sep='\t', index=False, header=True)


@click.option(
    '--input-dir',
    help='GCS path to the directory containing the SusiE results',
    default='gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/fine_mapping/susie_mcv/output_files',
)
@click.option('--celltypes', help='Comma-separated list of cell types to concatenate')
@click.option('--pip-prune-threshold', type=float, default=0.5, help='PIP threshold for pruning results')
@click.command()
def main(input_dir, celltypes, pip_prune_threshold):
    b = get_batch(name='Concatenate associatr-atac results')
    for cell_type in celltypes.split(','):
        j = b.new_python_job(
            name=f'Concatenate associatr-atac results for {cell_type}',
        )
        j.cpu(1)
        j.storage('5G')

        j.call(
            concatenator,
            input_dir,
            cell_type,
            pip_prune_threshold
        )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
