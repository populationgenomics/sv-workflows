#!/usr/bin/env python3

"""
This script prepares eTR inputs from every cell type for mash (including the null model).Prerequisite to the concatenate_inputs.py script.

analysis-runner --dataset "tenk10k" \
    --description "Prepare inputs for mashr" \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --access-level "test" \
    --output-dir "str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed" \
    process_inputs.py
"""


import click
import pandas as pd

from cpg_utils.hail_batch import get_batch, output_path
from cpg_utils import to_path


def cell_chrom_parser(cell, chrom, estrs_coord_chrom, meta_input_dir):
    from cpg_utils.hail_batch import output_path
    import pandas as pd

    cell_df = pd.DataFrame()
    for index, row in estrs_coord_chrom.iterrows():
        chrom = row['chr']
        gene_name = row['gene_name']
        pos = row['pos']
        motif = row['motif']
        try:
            df = pd.read_csv(
                f'gs://{meta_input_dir}/{cell}/{chrom}/{gene_name}_100000bp_meta_results.tsv',
                sep='\t',
            )

            df = df[(df['pos'] == pos) & (df['motif'] == motif)]
            beta = df['coeff_meta_fixed'].iloc[0]
            se = df['coeff_meta_fixed'].iloc[0]
            # save the beta and se into a row:
            beta_se = pd.DataFrame(
                {
                    'chrom': chrom,
                    'pos': pos,
                    'motif': motif,
                    'gene': gene_name,
                    f'{cell}_beta': beta,
                    f'{cell}_se': se,
                },
                index=[0],
            )
            cell_df = pd.concat([cell_df, beta_se], axis=0)
        except:
            continue

    cell_df.to_csv(
        output_path(f'mashr/beta_se/{cell}/{chrom}/beta_se.tsv'),
        sep='\t',
        index=False,
    )


def cell_chrom_parser_null(cell, chrom, meta_input_dir):
    from cpg_utils.hail_batch import output_path
    import pandas as pd

    gene_files = list(to_path(f'{meta_input_dir}/{cell}/chr{chrom}').rglob('*.tsv'))
    master_df = pd.DataFrame()
    for gene_file in gene_files:
        gene_name = str(gene_file).split('/')[-1].split('_')[0]
        df = pd.read_csv(gene_file, sep='\t')
        df['gene'] = gene_name
        df[f'{cell}_beta'] = df['coeff_meta_fixed']
        df[f'{cell}_se'] = df['se_meta_fixed']
        df = df[['chr', 'pos', 'motif', 'ref_len', 'gene', f'{cell}_beta', f'{cell}_se']]
        master_df = pd.concat([master_df, df], axis=0)

    master_df.to_csv(
        output_path(f'mashr/chr2_null_beta_se/{cell}/chr{chrom}/beta_se.tsv'),
        sep='\t',
        index=False,
    )


@click.option(
    '--cell-types',
    default='CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC',
    help='File containing eQTLs passing FDR threshold',
)
@click.option(
    '--estrs-coord-path',
    default='gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/cell-type-spec/estrs.csv',
    help='File containing eSTRs coordinates',
)
@click.option(
    '--meta-input-dir',
    default='gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/meta_results/meta_with_fixed/v2',
)
@click.command()
def main(cell_types, estrs_coord_path, meta_input_dir):
    b = get_batch(name='Process inputs for mashr')
    celltypes = cell_types.split(',')
    # load in the list of esnps passing FDR <5% across all cell types:
    estrs_coord = pd.read_csv(estrs_coord_path)

    for cell in celltypes:
        for chrom in range(1, 23):
            estrs_coord_chrom = estrs_coord[estrs_coord['chr'] == f'chr{chrom}']
            if to_path(output_path(f'mashr/beta_se/{cell}/{chrom}/beta_se.tsv')).exists():
                continue
            job = b.new_python_job(f'Prep eTRs for mashr {cell} {chrom}')
            job.cpu(0.25)
            job.call(cell_chrom_parser, cell, chrom, estrs_coord_chrom, meta_input_dir)

        for chrom in [2]:
            if to_path(output_path(f'mashr/chr3_null_beta_se/{cell}/chr{chrom}/beta_se.tsv')).exists():
                continue
            job = b.new_python_job(f'Prep eTRs for mashr {cell} {chrom} null')
            job.cpu(0.25)
            job.call(cell_chrom_parser_null, cell, chrom, meta_input_dir)

    b.run(wait=False)


if __name__ == '__main__':
    main()
