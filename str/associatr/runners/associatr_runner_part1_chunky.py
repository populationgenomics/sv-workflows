#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script is step 1 of 3 for running associaTR.
It aims to:
 - intersect genes found in the scRNA dataset (cell type + chr specific) with GENCODE v42 annotations.
 - write out pseudobulk and covariates files for each cell type

 Use this for larger cell files as it subsets for the genes of interest first before pseudobulk aggregation.

 analysis-runner --dataset "tob-wgs" \
    --description "prepare expression files for associatr" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
     associatr_runner_part1_chunky.py  --celltypes=CD4_NC_part2 --chromosomes=chr22

"""
import json
import click
import pandas as pd
import hail as hl
import hailtop.batch as hb
from cpg_utils import to_path
import numpy as np
import logging
logging.basicConfig(level = logging.INFO)
import csv
import sys

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

def build_pseudobulk(celltype, chromosomes):
    #intersect genes with gencode v42
    gencode = pd.read_table("gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/gencode.v42.annotation.gff3.gz", comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])

    gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)
    gencode_genes["gene_name"],gencode_genes["ENSG"], = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))
    print(f'gencode_genes df memory: %s', gencode_genes.memory_usage(deep=True).sum())


    #subset gencode annotation file for relevant chromosome(s)
    chromosome_list = chromosomes.split(',')
    gencode_genes = gencode_genes[gencode_genes['seqname'].isin(chromosome_list)]


    #select gene names in the sc-input file that match gencode_genes
    csv.field_size_limit(sys.maxsize)
    with to_path(f'gs://cpg-tob-wgs-test/hoptan-str/associatr/sc-input/{celltype}_sc_pheno_cov.tsv').open('r') as f:
        reader = csv.reader(f)
        header = next(reader)  # get the header row
    header = header[0].split('\t')
    genes_in_chr = list(set(header).intersection(set(gencode_genes['gene_name'])))
    genes_in_chr+=['individual','sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2']


    print(f'starting..')
    pheno_cov_file_path = f'gs://cpg-tob-wgs-test/hoptan-str/associatr/sc-input/{celltype}_sc_pheno_cov.tsv'
    pheno_cov_sc_input = pd.read_csv(pheno_cov_file_path, sep='\t', usecols=genes_in_chr)
    print(f'loaded in pheno_cov file for {celltype}')
    print('f pheno_cov_sc_input df memory: %s', pheno_cov_sc_input.memory_usage(deep=True).sum())
    covariates = pheno_cov_sc_input[['individual','sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2']]
    covariates = covariates.dropna() #removes rows with NaN values (likely a data processing error)
    pheno_sc_input = pheno_cov_sc_input.drop(columns=['barcode','sex','pc1','pc2','pc3','pc4','pc5','pc6','age','pf1','pf2'])
    sample_mapping_df = pd.read_csv('gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv', delimiter='\t')


    ## pseudo-bulk (mean gene counts, grouped by individual)
    # Group by 'individuals' and 'genes', then calculate the mean
    pseudobulk = pheno_sc_input.groupby('individual').mean()
    pseudobulk.reset_index(inplace=True)
    print(f'pseudobulk df memory: %s', pseudobulk.memory_usage(deep=True).sum())

    # Retain genes that are expressed in at least 10% of the individuals
    threshold = len(pseudobulk) * 0.1
    filtered_columns = [col for col in pseudobulk.columns if col != 'individual' and (pseudobulk[col] > 0).sum() >= threshold]
    filtered_columns+= ['individual']
    pseudobulk = pseudobulk[filtered_columns]
    print(f'pseudobulk df memory after filtering: %s', pseudobulk.memory_usage(deep=True).sum())

    #apply ln(x+1) transformation to the pseudobulk mean gene counts of the remaining genes
    individual_ids = pseudobulk['individual']
    pseudobulk_log_transformed = pseudobulk[pseudobulk.columns.difference(['individual'])].apply(lambda x: np.log1p(x))
    pseudobulk = pd.concat([individual_ids,pseudobulk_log_transformed], axis=1)
    print(f'pseudo-bulk df memory after log transformation: %s', pseudobulk.memory_usage(deep=True).sum())

    #add CPG IDs to pseudobulk dataframe
    pseudobulk = pseudobulk.merge(sample_mapping_df, left_on='individual', right_on = "OneK1K_ID", how='inner')
    pseudobulk['InternalID']=pseudobulk['InternalID'].str[3:] #removes 'CPG' prefix because associaTR requires strictly numeric IDs
    pseudobulk['InternalID'] = pseudobulk['InternalID'].astype(float)

    ##write pseudobulk to GCS
    pseudobulk.to_csv(f'gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/pseudobulk/{celltype}_{chromosomes}_pseudobulk.tsv', sep='\t', index=False)

    ##write covariates to GCS
    covariates.to_csv(f'gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/covariates/{celltype}_covariates.tsv', sep='\t', index=False)

    #intersect genes with gencode v42
    for chromosome in chromosomes.split(','):
        #subset gencode annotation file for relevant chromosome
        gencode_genes_chr = gencode_genes[gencode_genes['seqname']==chromosome]

        # Extract the gene names from the 'gene_name' column of 'gencode' DataFrame
        gencode_gene_names = set(gencode_genes_chr['gene_name'])

        # Get the genes with scRNA data
        gene_columns = [col for col in pseudobulk.columns if col != 'individual']

        # Create a set of gene names (representing gens with scRNA data), by iterating through columns
        pseudobulk_gene_names = set()
        for col in gene_columns:
            pseudobulk_gene_names.add(col)

        # Find the intersection of gene names between genes with scRNA data and GENCODE gene names
        pseudobulk_gene_names = gencode_gene_names.intersection(pseudobulk_gene_names)

        # write genes array to GCS directly
        with to_path(f'gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/scRNA_gene_lists/{celltype}/{chromosome}_{celltype}_filtered_genes.json').open('w') as write_handle:
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
        pseudobulk_job.memory('64G')
        pseudobulk_job.storage('20G')
        pseudobulk_job.cpu(4)
        pseudobulk_job.image(config['workflow']['driver_image'])
        pseudobulk_job.call(build_pseudobulk,celltype, chromosomes)

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter







