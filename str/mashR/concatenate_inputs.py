#!/usr/bin/env python3

"""
This script concatenates outputs from `process_inputs.py` for input into mashR to produce output files for chr1-22 per cell type.

analysis-runner --dataset "tenk10k" \
    --description "Prepare inputs for mashr" \
    --memory 64G \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --access-level "test" \
    --output-dir "potato" \
    concatenate_inputs.py



"""


import pandas as pd
import click

from cpg_utils.hail_batch import get_batch
from cpg_utils import to_path


@click.option(
    '--cell-types',
    default='CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC',
    help='File containing eQTLs passing FDR threshold',
)
@click.option(
    '--mash-process-inputs-dir',
    default='gs://cpg-tenk10k-test/str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed/mashr',
)
@click.command()
def main(cell_types, mash_process_inputs_dir):

    b = get_batch(name='Concatenate files for mashR')

    celltypes = cell_types.split(',')

    ## Concatenate beta and se files for each cell type across chromosomes
    for cell in celltypes:
        master_df = pd.DataFrame()
        for chrom in range(1, 23):
            df = pd.read_csv(
                f'{mash_process_inputs_dir}/beta_se/{cell}/chr{chrom}/beta_se.tsv',
                sep='\t',
            )
            if master_df.empty:
                master_df = df
            else:
                master_df = pd.concat([master_df, df], axis=0)
        master_df.to_csv(
            f'{mash_process_inputs_dir}/beta_se/all_chrom/{cell}_beta_se.tsv',
            sep='\t',
            index=False,
        )

    ## Concatenate beta and se files for all cell types across all chromosomes
    master_df = pd.DataFrame()
    for cell in celltypes:
        df = pd.read_csv(f'{mash_process_inputs_dir}/beta_se/all_chrom/{cell}_beta_se.tsv', sep='\t')

        if master_df.empty:
            master_df = df
        else:
            master_df = master_df.merge(df, on=['chrom', 'pos', 'motif', 'gene'], how='inner') # variants need to be shared across all cell types
    master_df.to_csv(
        f'{mash_process_inputs_dir}/beta_se/all_chrom/all_celltypes_beta_se.tsv',
        sep='\t',
        index=False,
    )


    b.run(wait=False)


if __name__ == '__main__':
    main()
