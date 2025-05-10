#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
import click

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

"""

analysis-runner --dataset "bioheart" --description "pval<0.05 filter" --access-level "test" \
    --output-dir "str/associatr/trans_pilot/tob_n950_bioheart_n975" concat.py --input-dir=gs://cpg-bioheart-test-analysis/str/associatr/trans_pilot/tob_n950_bioheart_n975/filtered_meta_results \
    --cell-types=CD4_TCM --chromosomes=22

    1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21


"""


def concat_and_filter(gene_files, cell_type, chromosome):
    """
    Concatenate all gene_files into one dataframe.
    Write one TSV per cell_type/chromosome pair.
    """
    import pandas as pd
    from cpg_utils.hail_batch import output_path
    all_dfs = []
    for gene_file in gene_files:
        df = pd.read_csv(gene_file, sep='\t')
        df['gene'] = gene_file.split('/')[-1].split('.')[0]
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)

    gcs_output = output_path(f'filtered_meta_results/{cell_type}/chr{chromosome}/combined_filtered.tsv', 'analysis')
    with to_path(gcs_output).open('w') as f:
        combined_df.to_csv(f, sep='\t', index=False)


@click.option(
    '--input-dir',
    help='GCS path to the raw meta results of associaTR',
    type=str,
)
@click.option(
    '--cell-types',
    help='Comma-separated list of cell types',
)
@click.option(
    '--chromosomes',
    help='Comma-separated list of chromosomes',
)
@click.command()
def main(input_dir, cell_types, chromosomes):
    """
    Concatenates and filters raw meta-analysis results for multiple cell types and chromosomes.
    """
    for cell_type in cell_types.split(','):
        b = get_batch(name=f'Concat {cell_type}')
        for chromosome in chromosomes.split(','):
            gene_files = list(to_path(f'{input_dir}/{cell_type}/chr{chromosome}').glob('*.tsv'))

            j = b.new_python_job(
                name=f'Concat for {cell_type} chr{chromosome}',
            )

            j.call(concat_and_filter, gene_files, cell_type, chromosome)
    b.run(wait=False)


if __name__ == '__main__':
    main()
