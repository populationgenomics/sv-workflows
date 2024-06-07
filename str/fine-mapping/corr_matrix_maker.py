#!/usr/bin/env python3

"""
This file is used to create a correlation matrix between STR and SNP genotypes, which is required by fine-mapping methods.

Output: Correlation matrix, and list of variants (as appears in the matrix)

Workflow (for each cell type):
1) Extract genes where the eSTR is significant (default = FDR<0.05).
For each gene,
2) Extract the STR and SNP coordinates where the association signal is p < 4e-4, to reduce computational burden of fine-mapper.
3) Obtain the genotypes for each extracted STR and SNP in 2)
4) Calculate the correlation matrix between STR and SNP genotypes.

association-runner --dataset "bioheart" \
    --description "Calculate LD between STR and SNPs" \
    --access-level "test" \
    --output-dir "str/associatr/fine_mapping/prep_files" \
    ld_runner.py --snp-vcf-dir=gs://cpg-bioheart-test/str/dummy_snp_vcf\
    --str-vcf-dir=gs://cpg-bioheart-test/str/saige-qtl/input_files/vcf/v1-chr-specific \
    --celltypes=ASDC \
    --associatr-dir=gs://cpg-bioheart-test/str/associatr/snps_and_strs/tob_n1055_and_bioheart_n990/meta_results


"""

import ast
import numpy as np
import click
import pandas as pd

from hailtop.batch import ResourceGroup

from cpg_utils import to_path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch

def ld_parser(
    snp_vcf_path: ResourceGroup,
    str_vcf_path: ResourceGroup,
    associatr_results: pd.DataFrame,
    gene: str,
    celltype: str,
) -> str:
    import pandas as pd
    from cyvcf2 import VCF

    # create empty DF to store the relevant GTs (SNPs, STRs)
    snp_df = pd.DataFrame(columns=['individual'])
    str_df = pd.DataFrame(columns=['individual'])

    # cyVCF2 reads the SNP,STR VCF
    snp_vcf = VCF(snp_vcf_path['vcf'])
    snp_df['individual'] = snp_vcf.samples
    str_vcf = VCF(str_vcf_path['vcf'])
    str_df['individual'] = str_vcf.samples

    for index, row in associatr_results.iterrows():
        pos = row['pos']
        motif = row['motif']
        if '-' in row['motif']: # SNP
            chr_num = row['chr'][3:]
            for variant in snp_vcf(f'{chr_num}:{pos}-{pos}'):
                gt = variant.gt_types
                gt[gt == 3] = 2
                snp_vcf[f'chr{chr_num}:{pos}_{motif}'] = gt
                break
        else: # STR
            chrom = row['chr']
            for variant in str_vcf(f'{chrom}:{pos}-{pos}'):
                if str(variant.INFO.get('RU')) == motif:
                    genotypes = variant.format('REPCN')
                    # Replace '.' with '-99/-99' to handle missing values
                    genotypes = np.where(genotypes == '.', '-99/-99', genotypes)

                    # Split each element by '/'
                    split_genotypes = [genotype.split('/') for genotype in genotypes]

                    # Convert split_genotypes into a numpy array for easier manipulation
                    split_genotypes_array = np.array(split_genotypes)

                    # Convert the strings to integers and sum them row-wise
                    sums = np.sum(split_genotypes_array.astype(int), axis=1)
                    # set dummy -198 value to np.nan
                    sums = np.where(sums == -198, np.nan, sums)

                    str_df[f'chr{chr_num}:{pos}_{motif}'] = sums
                    break
    # merge the two dataframes
    merged_df = str_vcf.merge(snp_vcf, on='individual')

    # calculate pairwise correlation of every variant
    corr_matrix = merged_df.drop(columns='individual').corr()
    corr_matrix.to_csv(f'correlation_matrix/{celltype}/{gene}_correlation_matrix.tsv', sep='\t')
    corr_matrix.columns.to_csv(f'correlation_matrix/{celltype}/{gene}_correlation_matrix_variants.tsv', sep='\t')

    print("Wrote correlation matrix to bucket")

@click.option(
    '--snp-vcf-dir',
    help='GCS file dir to SNP VCF files.',
    type=str,
)
@click.option(
    '--str-vcf-dir',
    help='GCS file dir to STR VCF files.',
    type=str,
)
@click.option(
    '--fdr-cutoff',
    help='FDR cutoff to use for eGene selection',
    type=float,
    default=0.05,)
@click.option(
    '--pval-cutoff',
    help='P-value cutoff to use for associatr results (reduce finemapping computational burden)',
    default=4e-4,
)
@click.option('--celltypes', help='Cell types to use for coloc', type=str)

@click.option(
    '--str-fdr-dir',
    help='Path to STR FDR dir',
    default='gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat',
)

@click.option('--associatr-dir', help='Path to STR-SNP associatr results', default = 'gs://cpg-bioheart-main-analysis/str/associatr/snps_and_strs/tob_n1055_and_bioheart_n990/meta_results')

@click.option('--job-cpu', default=1)
@click.option('--job-storage', default='20G')
@click.command()
def main(
    snp_vcf_dir: str,
    str_vcf_dir: str,
    celltypes: str,
    fdr_cutoff: float,
    str_fdr_dir: str,
    job_cpu: int,
    job_storage: str,
    associatr_dir: str,
    pval_cutoff: float,
):
    b = get_batch()
    for celltype in celltypes.split(','):
        # read in STR eGene annotation file
        str_fdr_file = f'{str_fdr_dir}/{celltype}_qval.tsv'
        str_fdr = pd.read_csv(str_fdr_file, sep='\t')
        str_fdr = str_fdr[str_fdr['qval'] < fdr_cutoff]  # subset to eGenes passing FDR 5% threshold by default
        for index, row in str_fdr.iterrows():
            gene = row['gene_name']
            chrom = ast.literal_eval(row['chr'].iloc[0])
            #test only
            if chrom!= 'chr20':
                continue
            try:
                associatr_file = f'{associatr_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv'
                associatr = pd.read_csv(associatr_file, sep='\t')
            except:
                f'No associatr results for this gene: {gene}'
                continue
            associatr = associatr[associatr['pval_meta'] < pval_cutoff]
            if associatr.empty:
                f'No associatr results for this gene: {gene}'
                continue
            snp_vcf_path = f'{snp_vcf_dir}/{chrom}_common_variants.vcf.bgz'
            str_vcf_path = f'{str_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'
            # run coloc
            ld_job = b.new_python_job(
                f'LD calc for {gene}: {celltype}',
            )
            ld_job.cpu(job_cpu)
            ld_job.storage(job_storage)
            snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_path, 'csi': snp_vcf_path + '.csi'})
            str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'csi': str_vcf_path + '.csi'})

            ld_job.call(
                ld_parser,
                snp_input,
                str_input,
                associatr,
                gene,
                celltype,
            )

    b.run(wait=False)
if __name__ == '__main__':
    main()