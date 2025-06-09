#!/usr/bin/env python3

"""

This script is the first step in assessing cell type-specificity of eQTLs.
It prepares input files for the next step, which is running the meta-analysis; and also outputs a file containing eQTLs with opposite signed betas for each cell type.

analysis-runner --dataset "tenk10k" --description "eqtl_file_prep" --access-level "test" \
--output-dir "str/associatr/cell-type-spec" file_prep.py --eqtl-file=gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/cell-type-spec/estrs.csv \
--associatr-dir=cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed/v2




"""

import logging

import click
import pandas as pd

from cpg_utils.hail_batch import get_batch


def meta_eqt_file_prep(cell_type_eqtls, cell_type, associatr_dir):
    import pandas as pd

    from cpg_utils.hail_batch import output_path

    meta_input_df = pd.DataFrame()
    opposite_signed_betas = pd.DataFrame()
    cell_type_list = 'CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,CD4_TCM_permuted,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC'
    cell_type_array = cell_type_list.split(',')

    for index, row in cell_type_eqtls.iterrows():
        gene = row['gene_name']
        chrom = row['chr']
        pos = row['pos']
        end = row['end']
        motif = row['motif']
        pval = row['pval_meta_fixed']
        for cell_type2 in cell_type_array:
            if cell_type2 != cell_type:
                file = f'{associatr_dir}/{cell_type2}/{chrom}/{gene}_100000bp_meta_results.tsv'
                try:
                    eqtl_df2 = pd.read_csv(file, sep='\t')
                    eqtl_df2 = eqtl_df2[eqtl_df2['pos'] == pos]
                    eqtl_df2['motif_len'] = eqtl_df2['motif'].str.len()
                    eqtl_df2['end'] = (
                        (
                            eqtl_df2['pos'].astype(float)
                            + eqtl_df2['ref_len'].astype(float) * eqtl_df2['motif_len'].astype(float)
                        )
                        .round()
                        .astype(int)
                    )
                    eqtl_df2 = eqtl_df2[eqtl_df2['end'] == end]
                    eqtl_df2 = eqtl_df2[eqtl_df2['motif'] == motif]
                    eqtl_df2_coeff = eqtl_df2['coeff_meta_fixed'].iloc[0]
                    eqtl_df2_se = eqtl_df2['se_meta_fixed'].iloc[0]

                    ## add a row to the meta_input_df
                    new_row = pd.DataFrame(
                        {
                            'chrom': [chrom],
                            'pos': [pos],
                            'end': [end],
                            'motif': [motif],
                            'gene_name': [gene],
                            'celltype_main': [cell_type],
                            'coeff_main': [row['coeff_meta_fixed']],
                            'se_main': [row['se_meta_fixed']],
                            'pval_main': [pval],
                            'cell_type2': [cell_type2],
                            'coeff_2': [eqtl_df2_coeff],
                            'se_2': [eqtl_df2_se],
                        },
                    )
                    meta_input_df = pd.concat([meta_input_df, new_row], ignore_index=True)
                    if row['coeff_meta_fixed'] * eqtl_df2_coeff < 0:
                        opposite_signed_betas = pd.concat([opposite_signed_betas, new_row], ignore_index=True)
                except FileNotFoundError:
                    logging.info(f'File {file} not found')
                    continue

    o_file_path = output_path(f'prep_files/{cell_type}/meta_input_df.csv')
    o_file_path_opposite = output_path(f'prep_files/{cell_type}/opposite_signed_betas.csv')
    meta_input_df.to_csv(o_file_path, index=False)
    opposite_signed_betas.to_csv(o_file_path_opposite, index=False)


@click.option('--eqtl-file', help='File containing eQTLs passing FDR threshold')
@click.option('--associatr-dir', help='Directory containing associaTR raw outputs')
@click.command()
def main(eqtl_file, associatr_dir):
    df = pd.read_csv(eqtl_file)
    for cell_type in df['cell_type'].unique():
        cell_type_eqtls = df[df['cell_type'] == cell_type]
        j = get_batch(name='meta_eqt_file_prep').new_python_job(name=f'{cell_type}_meta_eqt_file_prep')
        j.call(meta_eqt_file_prep, cell_type_eqtls, cell_type, associatr_dir)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
