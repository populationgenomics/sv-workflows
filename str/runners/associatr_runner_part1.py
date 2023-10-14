#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script is step 1 of 2 for running associaTR.
It aims to:
 - intersect genes found in the scRNA dataset (cell type + chr specific) with GENCODE v42 annotations.
 - output cis window files for each gene in the above intersection

 analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files for associatr" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
     associatr_runner_part1.py  --celltypes=Plasma --chromosomes=chr22

"""
import json
import click
import pandas as pd
import hail as hl
import hailtop.batch as hb
from cpg_utils import to_path

from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch

from cpg_utils.hail_batch import output_path, init_batch, dataset_path

config = get_config()

def gene_info(x):
    # Extract ENSG and gene_level of evidence
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    g_id = list(filter(lambda x: 'gene_id' in x,  x.split(";")))[0].split("=")[1]
    g_id = g_id.split('.')[0] #removes the version number from ENSG ids
    return (g_name,g_id)

def pseudobulk(celltype, chromosomes):
    pheno_cov_file_path = f'gs://cpg-tob-wgs-test/hoptan-str/associatr/sc-input/{celltype}_sc_pheno_cov.tsv'
    pheno_cov_sc_input = pd.read_csv(pheno_cov_file_path, sep='\t')
    pheno_sc_input = pheno_cov_sc_input.drop(columns=['barcode','sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2'])

    ## pseudo-bulk (mean gene counts, grouped by individual)
    # Group by 'individuals' and 'genes', then calculate the mean
    pseudobulk = pheno_sc_input.groupby('individual').mean()
    pseudobulk.reset_index(inplace=True)

    # Basic filter - remove genes where all entries are either 0 or NA in the pseudobulk data frame
    mask = (pseudobulk == 0) | pseudobulk.isna()
    mask = mask.all(axis=0)
    pseudobulk = pseudobulk.loc[:, ~mask]

    #intersect genes with gencode v42
    gencode = pd.read_table("gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/gencode.v42.annotation.gff3.gz", comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
    gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)
    gencode_genes["gene_name"],gencode_genes["ENSG"], = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))

    for chromosome in chromosomes.split(','):
        #subset gencode annotation file for relevant chromosome
        gencode_genes = gencode_genes[gencode_genes['seqname']==chromosome]

        # Extract the gene names from the 'gene_name' column of 'gencode' DataFrame
        gencode_gene_names = set(gencode_genes['gene_name'])

        # Get the genes with scRNA data
        gene_columns = [col for col in pseudobulk.columns if col != 'individual']

        # Create a set of gene names (representing gens with scRNA data), by iterating through columns
        pseudobulk_gene_names = set()
        for col in gene_columns:
            pseudobulk_gene_names.add(col)

        # Find the intersection of gene names between genes with scRNA data and GENCODE gene names
        pseudobulk_gene_names = gencode_gene_names.intersection(pseudobulk_gene_names)

        # write genes array to GCS directly
        with to_path(output_path(f'input_files/scRNA_gene_lists/{celltype}/{chromosome}_{celltype}_filtered_genes.json')).open('w') as write_handle:
            json.dump(list(pseudobulk_gene_names), write_handle)

# inputs:
@click.option('--celltypes')
@click.option('--chromosomes', help=' eg chr22')
@click.command()
def main(
    celltypes, chromosomes
):
    """
    Run associatr_runner_part1.py
    """
    config = get_config()
    b = get_batch()
    init_batch()

    for celltype in celltypes.split(','):
        pseudobulk_job = b.new_python_job(name=f'Build pseudobulk and filter for {celltype}')
        pseudobulk_job.image(config['workflow']['driver_image'])
        pseudobulk_job.call(pseudobulk,celltype, chromosomes)

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter







