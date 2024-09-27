#!/usr/bin/env python3

"""

This script is the first step in assessing cell type-specificity of eQTLs.
It performs pairwise meta-analysis of eQTLs (e.g. testing maximum of two cell types) for each eQTL passing MTC (eg FDR <5%).

analysis-runner --dataset "bioheart" --description "eqtl_file_prep" --access-level "test" \
--output-dir "str/associatr/cell-type-spec" file_prep.py --eqtl-file=gs://cpg-bioheart-test/str/associatr/cell-type-spec/estrs.csv \
--associatr-dir=gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results




"""

import click
import pandas as pd

from cpg_utils.hail_batch import get_batch, output_path


def meta_eqt_file_prep(cell_type_eqtls, cell_type, associatr_dir):
    meta_input_df = pd.DataFrame()
    cell_type_list = 'CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,CD4_TCM_permuted,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC'
    cell_type_array = cell_type_list.split(',')

    for row in cell_type_eqtls.iterrows():
        gene = row['gene_name']
        chrom = row['chr']
        pos = row['pos']
        end = row['end']
        motif = row['motif']
        pval = row['pval']
        for cell_type2 in cell_type_array:
            if cell_type2 != cell_type:
                file = f'{associatr_dir}/{cell_type}/{chrom}/{gene}_100000bp_meta_results.tsv'
                try:
                    eqtl_df2 = pd.read_csv(file, sep='\t')
                except FileNotFoundError:
                    continue
                eqtl_df2 = eqtl_df2[eqtl_df2['gene'] == gene]
                eqtl_df2 = eqtl_df2[eqtl_df2['chr'] == chrom]
                eqtl_df2 = eqtl_df2[eqtl_df2['pos'] == pos]
                eqtl_df2 = eqtl_df2[eqtl_df2['end'] == end]
                eqtl_df2 = eqtl_df2[eqtl_df2['motif'] == motif]

                ## add a row to the meta_input_df
                new_row = pd.DataFrame(
                    {
                        'chrom': [chrom],
                        'pos': [pos],
                        'end': [end],
                        'motif': [motif],
                        'celltype_main': [cell_type],
                        'coeff_main': [row['coeff']],
                        'se_main': [row['se']],
                        'pval_main': [pval],
                        'cell_type2': [cell_type2],
                        'coeff_2': [eqtl_df2['pval']],
                        'se_2': [eqtl_df2['se']],
                    },
                )
                meta_input_df = pd.concat([meta_input_df, new_row], ignore_index=True)

    o_file_path = output_path(f'prep_files/{cell_type}/meta_input_df.csv')
    meta_input_df.to_csv(o_file_path, index=False)


@click.option('--eqtl-file', help='File containing eQTLs passing FDR threshold')
@click.option('--associatr-dir', help='Directory containing associaTR raw outputs')
@click.command()
def main(eqtl_file, associatr_dir):
    df = pd.read_csv(eqtl_file)
    df = df.drop_duplicates(subset=['chr', 'pos', 'motif'])
    #for cell_type in df['cell_type'].unique():
    for cell_type in ['ASDC']:
        cell_type_eqtls = df[df['cell_type'] == cell_type]
        j = get_batch(name='meta_eqt_file_prep').new_python_job(name=f'{cell_type}_meta_eqt_file_prep')
        j.call(meta_eqt_file_prep, cell_type_eqtls, cell_type, associatr_dir)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
