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
     associatr_runner_part1.py  --celltypes=Plasma --chromosomes=chr22 \

"""
import click
import pandas as pd
import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path

from cpg_utils.hail_batch import output_path

config = get_config()

def gene_info(x):
    # Extract ENSG and gene_level of evidence
    g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
    g_id = list(filter(lambda x: 'gene_id' in x,  x.split(";")))[0].split("=")[1]
    g_id = g_id.split('.')[0] #removes the version number from ENSG ids
    return (g_name,g_id)

def get_gene_cis_file(chromosome:str, gene: str, window_size: int, ofile_path: str):
    """Get gene cis window file"""
    gencode = pd.read_table("gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/gencode.v42.annotation.gff3.gz", comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
    gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)
    gencode_genes["gene_name"],gencode_genes["ENSG"] = zip(*gencode_genes.attribute.apply(gene_info))

    #subset gencode annotation file for relevant chromosome
    gencode_genes = gencode_genes[gencode_genes['seqname']==chromosome]

    #subsets gencode annotation file for relevant gene
    gene_info_gene = gencode_genes[gencode_genes['gene_name'] == gene]

    init_batch()
    # get chromosome, coordinates
    chrom = gene_info_gene['seqname'].values[0]
    start_coordinate = gene_info_gene['start'].values[0]
    end_coordinate = gene_info_gene['end'].values[0]
    # get gene body position (start and end) and add window
    left_boundary = max(1, start_coordinate - window_size)
    right_boundary = min(
        end_coordinate + window_size,
        hl.get_reference('GRCh38').lengths[chrom]
    )
    data = {'chromosome': chrom, 'start': left_boundary, 'end': right_boundary}
    pd.DataFrame(data, index=[gene]).to_csv(ofile_path, sep='\t', header=False)


# inputs:
@click.option('--celltypes')
@click.option('--chromosomes', help=' eg chr22')
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=50,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
@click.option('--cis-window-size', type=int,default=100000)
@click.command()
def main(
    celltypes, chromosomes, max_parallel_jobs, cis_window_size
):
    """
    Run associaTR processing pipeline
    """
    config = get_config()
    b = get_batch()
    init_batch()

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    for celltype in celltypes:
        pheno_cov_sc_input = pd.read_csv(f'gs://cpg-tob-wgs-test/sc-input/{celltype}_sc_pheno_cov.tsv', sep='\t')
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
        gencode_genes["gene_name"],gencode_genes["ENSG"], gencode_genes["gene_level"] = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))

        for chromosome in chromosomes:
            #subset gencode annotation file for relevant chromosome
            gencode_genes = gencode_genes[gencode_genes['seqname']=='chr22']

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

            for gene in pseudobulk_gene_names:
                # get gene cis-window file
                gene_cis_job = b.new_python_job(name=f'Build cis window files for {gene} [{celltype};{chromosome}]')
                gene_cis_job.image(config['workflow']['driver_image'])
                gene_cis_job.call(
                        get_gene_cis_file,chromosome,gene,cis_window_size,gene_cis_job.ofile
                    )
                b.write_output(gene_cis_job.ofile, output_path(f'input_files/cis_window_files/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.bed'))
                manage_concurrency_for_job(gene_cis_job)
    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter







