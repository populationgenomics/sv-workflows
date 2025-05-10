#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
This script filters meta_results for rows where P<0.05.

    analysis-runner --dataset "bioheart" --description "pval<0.05 filter" --access-level "test" \
    --output-dir "str/associatr/trans_pilot/tob_n950_bioheart_n975" nominal_pval_filter.py --input-dir=gs://cpg-bioheart-test-analysis/str/associatr/trans_pilot/tob_n950_bioheart_n975/meta_results \
    --cell-types=CD4_TCM --chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21

"""

import click

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path

def pval_filter(gene_file, cell_type, chromosome):
    """
    Filter the meta results for rows where P<0.05
    """
    import pandas as pd
    from cpg_utils.hail_batch import output_path

    # read the raw results
    gene_results = pd.read_csv(gene_file, sep='\t')
    # filter for pval < 0.05
    gene_results = gene_results[gene_results['pval_meta'] < 0.05]
    # write to output file
    gcs_output = output_path(f'filtered_meta_results/{cell_type}/chr{chromosome}/{gene_file.name}', 'analysis')
    with to_path(gcs_output).open('w') as f:
        gene_results.to_csv(f, sep='\t', index=False)

@click.option(
    '--input-dir',
    help='GCS path to the raw meta results of associaTR',
    type=str,
)
@click.option(
    '--cell-types',
    help='Name of the cell type, comma separated if multiple',
)
@click.option(
    '--chromosomes',
    help='Chromosome number eg 1, comma separated if multiple',
)
@click.command()
def main(input_dir, cell_types, chromosomes):
    """
    Extracts the raw p-values from the results of associaTR into one text file per cell type.
    """
    for cell_type in cell_types.split(','):
        b = get_batch(name = f'Pval filter for {cell_type}')
        for chromosome in chromosomes.split(','):
            gene_files = list(to_path(f'{input_dir}/{cell_type}/chr{chromosome}').glob('*.tsv'))
            for gene_file in gene_files:
                if to_path(output_path(f'filtered_meta_results/{cell_type}/chr{chromosome}/{gene_file.name}', 'analysis')).exists():
                    continue
                j = b.new_python_job(
                    name=f'Pval filter for {cell_type} {chromosome}: {gene_file.name}',
                )
                j.cpu(0.25).memory('lowmem')
                j.call(pval_filter, gene_file, cell_type, chromosome)
    b.run(wait=False)


if __name__ == '__main__':
    main()
