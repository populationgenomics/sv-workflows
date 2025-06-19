#!/usr/bin/env python3

import click

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path
import pandas as pd

"""
This script creates a BED file of all coordinates from a directory.

analysis-runner --dataset "tenk10k" --description "Create a dataframe from all cell types" --access-level "test" \
    --output-dir "str/associatr/fine_mapping/susie_finemap" \
    df_to_bed.py

"""


@click.option(
    '--input-dir',
    type=str,
    required=True,
    help='Directory containing the meta results',
    default='gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed/v2/min_p_0.05',
)
@click.command()
def main(input_dir):
    cell_types = 'B_intermediate,ILC,Plasmablast,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,NK,CD8_TEM,CD4_Naive,B_naive,CD4_TCM,gdT,dnT'
    b = get_batch()
    bed_entries = set()

    # Split the string into a list of characters
    cell_type_list = cell_types.split(',')
    for cell in cell_type_list:
        df = pd.read_csv(f'{input_dir}/{cell}/all_genes_pval_0.05.tsv', sep='\t')
        for row in df.itertuples(index=False):
            # Add as tuple to set
            bed_entries.add((row.chr, row.pos, row.end))

    # Convert set to sorted DataFrame
    bed_df = pd.DataFrame(list(bed_entries), columns=['chr', 'start', 'end'])
    bed_df = bed_df.sort_values(['chr', 'start', 'end'])

    # Ensure 'chr' prefix if not already
    bed_df['chr'] = bed_df['chr'].astype(str)
    bed_df['chr'] = bed_df['chr'].apply(lambda c: c if c.startswith('chr') else f'chr{c}')

    # Write to BED file
    bed_df.to_csv(f'{input_dir}/all_coordinates.bed', sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
