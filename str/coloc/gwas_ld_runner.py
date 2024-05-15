#!/usr/bin/env python3

"""

This script calculates pairwise LD (correlation coefficient) between the top eSTR and the lead SNP (per gene) in the GWAS.
(DOES NOT REQUIRE coloc_runner.py and coloc_results_parser.py to be run previously, nor does it run `coloc`.)

Workflow:
1) Extract genes that are associated with the eSTRs (FDR <5%)
2) Extract the SNP GWAS data for the cis-window (gene +/- 100kB)
3) Select the lead SNP from 2) (ie SNP with the lowest p-value)
4) Calculate pairwise correlation of the lead eSTR locus with the lead SNP.

analysis-runner --dataset "bioheart" \
    --description "Calculate LD between STR and SNPs" \
    --access-level "test" \
    --cpu=1 \
    --output-dir "str/associatr/freeze_1/gwas_ld/bioheart-only-snps" \
    gwas_ld_runner.py --snp-vcf-dir=gs://cpg-bioheart-test/str/dummy_snp_vcf \
    --str-vcf-dir=gs://cpg-bioheart-test/str/saige-qtl/input_files/vcf/v1-chr-specific \
    --gwas-file=gs://cpg-bioheart-test/str/gwas_catalog/hg38.EUR.IBD.gwas_info03_filtered.assoc_for_gwas_ld.csv \
    --celltypes=CD4_TCM \
    --phenotype=ibd

"""

import ast

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch


# Function to process each element in the 'chr' column
def process_chr_element(element):
    try:
        # Safely evaluate the string to a list
        chr_list = ast.literal_eval(element)
        # Ensure the list is not empty and the first element is a string
        if chr_list and isinstance(chr_list[0], str):
            # Slice the string starting from the fourth character
            return chr_list[0][3:]
    except (ValueError, SyntaxError):
        # Handle cases where ast.literal_eval fails
        return None


def ld_parser(
    snp_vcf_path,
    str_vcf_path,
    phenotype,
    celltype,
    gene_annotation_file,
    str_fdr_dir,
    gwas_file,
    chromosome,
):
    from cyvcf2 import VCF

    from cpg_utils import to_path

    # read in gwas catalog file
    gwas_catalog_orig = pd.read_csv(gwas_file)
    # read in gene_annotation_table
    gene_annotation_table = pd.read_csv(gene_annotation_file)
    # load in the str fdr file
    str_fdr_file = f'{str_fdr_dir}/{celltype}_qval.tsv'
    str_fdr = pd.read_csv(str_fdr_file, sep='\t')
    str_fdr = str_fdr[str_fdr['qval'] < 0.05]  # subset to eGenes passing FDR 5% threshold

    # Apply the function to the 'chr' column and create a new 'chrom' column
    str_fdr['chrom_num'] = str_fdr['chr'].apply(process_chr_element)
    str_fdr = str_fdr[str_fdr['chrom_num'] == str(chromosome)]  # subset to the chromosome

    for gene in str_fdr['gene_name']:
        print(f'Processing gene: {gene}')

        # obtain snp cis-window coordinates for the gene
        gene_table = gene_annotation_table[gene_annotation_table['gene_ids'] == gene]  # subset to particular ENSG ID
        start_snp_window = float(gene_table['start'].iloc[0]) - 100000  # +-100kB window around gene
        end_snp_window = float(gene_table['end'].iloc[0]) + 100000
        chrom = gene_table['chr'].iloc[0][3:]
        print('Obtained SNP window coordinates')

        # subset the gwas catalog to the snp_window
        gwas_catalog = gwas_catalog_orig[gwas_catalog_orig['CHR'] == int(chrom)]
        gwas_catalog = gwas_catalog[gwas_catalog['BP'] >= start_snp_window]
        gwas_catalog = gwas_catalog[gwas_catalog['BP'] <= end_snp_window]

        if gwas_catalog.empty:
            print('No SNP GWAS data for ' + gene + ' in the cis-window: skipping....')
            continue

        # obtain lead SNP (lowest p-value) in the snp_window of the GWAS
        lowest_p_row = gwas_catalog.loc[gwas_catalog['P'].idxmin()]
        lead_snp_chr = lowest_p_row['CHR']
        lead_snp_bp = lowest_p_row['BP']
        lead_snp_end = lowest_p_row['BP'] + 1
        lead_snp_locus = f'{lead_snp_chr}:{lead_snp_bp}-{lead_snp_end}'

        # obtain top STR locus for the gene (if multiple are tied - iterate over each)
        str_fdr_gene = str_fdr[str_fdr['gene_name'] == gene]
        for estr in zip(
            ast.literal_eval(str_fdr_gene['chr'].iloc[0]),
            ast.literal_eval(str_fdr_gene['pos'].iloc[0]),
        ):
            chr_num = estr[0][3:]
            pos = estr[1]
            end = str(int(pos) + 1)
            str_locus = f'{chr_num}:{pos}-{end}'
            write_path = output_path(
                f'gwas_ld/{phenotype}/{celltype}/{gene}_chr{chr_num}_{pos}_gwas_ld_results.csv',
                'analysis',
            )

            if to_path(write_path).exists():
                print(f'GWAS LD for {gene} and {str_locus} already exists. Skipping...')
                continue

            print(f'Running LD for {gene} and {str_locus}')

            # create empty DF to store the relevant GTs (SNPs)
            df = pd.DataFrame(columns=['individual'])

            # cyVCF2 reads the SNP VCF
            vcf = VCF(snp_vcf_path['vcf'])
            df['individual'] = vcf.samples
            print('Reading SNP VCF with VCF()')

            print(f'Starting to subset VCF for the lead SNP {lead_snp_locus}')
            for variant in vcf(lead_snp_locus):
                geno = variant.gt_types  # extracts GTs as a numpy array
                locus = variant.CHROM + ':' + str(variant.POS)
                df_to_append = pd.DataFrame(geno, columns=[locus])  # creates a temp df to store the GTs for one locus
                break  # take the first SNP locus

            if df_to_append.empty:
                print(f'No GTs for SNP {locus} in the VCF, skipping...')
                continue

            # concatenate results to the main df
            df = pd.concat([df, df_to_append], axis=1)
            print("Finished subsetting VCF for lead SNP")

            # extract GTs for the one STR
            str_vcf = VCF(str_vcf_path['vcf'])
            for variant in str_vcf(str_locus):
                print(f'Captured STR with POS:{variant.POS}')
                ds = variant.format('DS')
                ds_list = []
                for i in range(len(ds)):
                    ds_list.append(ds[i][0])
                target_data = {'individual': str_vcf.samples, str_locus: ds_list}
                target_df = pd.DataFrame(target_data)
                break  # take the first STR locus

            # merge the two dataframes
            merged_df = df.merge(target_df, on='individual')

            # calculate pairwise correlation of every SNP locus with target STR locus
            correlation_series = merged_df.drop(columns='individual').corrwith(merged_df[str_locus])

            correlation_df = pd.DataFrame(correlation_series, columns=['correlation'])
            correlation_df['locus'] = correlation_df.index

            # drop the STR locus from the list of SNPs (it will automatically have a correlation of 1)
            correlation_df = correlation_df[correlation_df['locus'] != str_locus]

            # add some attributes
            correlation_df['gene'] = gene
            correlation_df['str_locus'] = str_locus
            correlation_df['celltype'] = celltype

            # write to output_path
            correlation_df.to_csv(write_path, index=False)


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
@click.option('--phenotype', help='Phenotype to use for coloc', type=str)
@click.option('--celltypes', help='Cell types to use for coloc', type=str)
@click.option(
    '--gene-annotation-file',
    help='Path to gene annotation file',
    default='gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv',
)
@click.option(
    '--str-fdr-dir',
    help='Path to STR FDR dir',
    default='gs://cpg-bioheart-test/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat',
)
@click.option(
    '--gwas-file',
    help='Path to GWAS catalog (ensure only three columns CHR, BP, and P)',
)
@click.option('--job-cpu', default=1)
@click.option('--job-storage', default='20G')
@click.command()
def main(
    snp_vcf_dir: str,
    str_vcf_dir: str,
    phenotype: str,
    celltypes: str,
    gene_annotation_file: str,
    str_fdr_dir: str,
    job_cpu: int,
    job_storage: str,
    gwas_file: str,
):
    b = get_batch(name = 'GWAS LD runner')

    for celltype in celltypes.split(','):
        for chromosome in [20]:
            ld_job = b.new_python_job(
                f'LD calc for chr{chromosome}; {celltype}',
            )
            ld_job.cpu(job_cpu)
            ld_job.storage(job_storage)

            snp_vcf_path = f'{snp_vcf_dir}/chr{chromosome}_common_variants_renamed.vcf.bgz'
            str_vcf_path = f'{str_vcf_dir}/hail_filtered_chr{chromosome}.vcf.bgz'

            snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_path, 'csi': snp_vcf_path + '.csi'})
            str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'csi': str_vcf_path + '.csi'})

            ld_job.call(
                ld_parser,
                snp_input,
                str_input,
                phenotype,
                celltype,
                gene_annotation_file,
                str_fdr_dir,
                gwas_file,
                chromosome,
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()
