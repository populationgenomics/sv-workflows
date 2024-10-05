#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script writes out genes for every cell type if the gene has a lead signal that is not a SNP.
We test all genes tested, not just eGenes that pass an FDR.
Used to generate one stat in the paper.

 analysis-runner  --dataset "bioheart" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr" \
python3 lead_tr_of_all_genes_tested.py

"""
import json

import pandas as pd

import hail as hl
import hailtop.batch as hb

from cpg_utils.hail_batch import get_batch


def gene_with_lead_tr_parser(
    chromosome,
    cell_type,
):
    """ """
    import numpy as np
    import scanpy as sc
    from scipy.stats import norm

    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    genes_with_lead_tr = []
    chromosome = f'chr{chromosome}'
    gene_file = f'gs://cpg-bioheart-test/str/associatr/input_files/240_libraries_tenk10kp1_v2/scRNA_gene_lists/1_min_pct_cells_expressed/{cell_type}/{chromosome}_{cell_type}_gene_list.json'
    with to_path(gene_file).open() as file:
        genes = json.load(file)

    for gene in genes:
        try:
            eqtl_results = pd.read_csv(
                f'gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results/{cell_type}/{chromosome}/{gene}_100000bp_meta_results.tsv',
                sep='\t',
            )
        except FileNotFoundError:
            print(f'No eQTL results found for {gene}... skipping')
            continue
        # get row(s) with minimum p-value
        min_pval = eqtl_results['pval_meta'].min()
        smallest_pval_rows = eqtl_results[eqtl_results['pval_meta'] == min_pval]
        # check if all rows are SNPs:
        all_motif_dash = smallest_pval_rows['motif'].str.contains('-').all()
        if all_motif_dash:
            print(f'Lead signal(s) is a SNP for {gene}... skipping')
            continue
        else:
            genes_with_lead_tr.append(gene)

    with to_path(
        output_path(
            f'genes_with_lead_trs/{cell_type}/{chromosome}_{cell_type}_gene_list.json',
        ),
    ).open('w') as write_handle:
        json.dump(genes_with_lead_tr, write_handle)


def main():
    """
    Get all genes that have a lead signal that is not a SNP.
    """
    b = get_batch(name='Extract genes with a tr lead signal')

    celltypes = 'gdT,B_intermediate,ILC,Plasmablast,dnT,ASDC,cDC1,pDC,NK_CD56bright,MAIT,B_memory,CD4_CTL,CD4_Proliferating,CD8_Proliferating,HSPC,NK_Proliferating,cDC2,CD16_Mono,Treg,CD14_Mono,CD8_TCM,CD4_TEM,CD8_Naive,CD4_TCM,NK,CD8_TEM,CD4_Naive,B_naive'
    celltypes = celltypes.split(',')
    for cell_type in ['ASDC']:
        # for chrom in range(1, 23):
        for chrom in [1]:
            j = b.new_python_job(
                name=f'Get genes with lead tr signal {cell_type}: {chrom}',
            )
            j.call(
                gene_with_lead_tr_parser,
                chrom,
                cell_type,
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments
