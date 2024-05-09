#!/usr/bin/env python3

"""

Assumes coloc_runner.py and coloc_results_parser.py have been run previously.

This script calculates pairwise LD (correlation coefficient) between one STR locus and all SNPs in a colocalised locus.

Workflow:
1)Extract genes that have >50% posterior probability of a shared causal variant from coloc results.
2) Extract the coordinates of the cis-window (gene +/- 100kB) for each gene using a gene annotation file.
3) Extract the top STR locus (ie passed FDR 5% threshold) for each gene in 1). May be multiple STRs per gene (if tied for smallest ACAT-corrected p-value).
4) Run LD parser() which:
- Extracts GTs for all SNPs (from a chr-specific VCF) in the cis-window for a gene.
- Extracts GTs for the specified STR locus associated with the gene.
- Calculates pairwise correlation of every SNP locus with the target STR locus.
- Save the SNP with the highest absolute correlation to a TSV file. Output to GCP

analysis-runner --dataset "bioheart" \
    --description "Calculate LD between STR and SNPs" \
    --access-level "full" \
    --cpu=1 \
    --output-dir "str/associatr/freeze_1/coloc_ld/bioheart-only-snps" \
    ld_runner.py --snp-vcf-dir=gs://cpg-bioheart-main/saige-qtl/bioheart_n990/input_files/genotypes/vds-bioheart1-0 \
    --str-vcf-dir=gs://cpg-bioheart-test/str/saige-qtl/input_files/vcf/v1-chr-specific \
    --coloc-dir=gs://cpg-bioheart-test/str/associatr/coloc \
    --phenotype=ibd \
    --celltypes=CD4_TCM

"""

import ast

import click
import pandas as pd

from hailtop.batch import ResourceGroup

from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch


def ld_parser(
    snp_vcf_path: ResourceGroup,
    str_vcf_path: ResourceGroup,
    str_locus: str,
    window: str,
    gwas_snp_path: str,
    gene: str,
    celltype: str,
) -> str:
    import pandas as pd
    from cyvcf2 import VCF

    # create empty DF to store the relevant GTs (SNPs)
    df = pd.DataFrame(columns=['individual'])

    # cyVCF2 reads the SNP VCF
    vcf = VCF(snp_vcf_path['vcf'])
    df['individual'] = vcf.samples
    print('Reading SNP VCF with VCF()')

    print('Starting to subset VCF for window...')
    for variant in vcf(window):
        geno = variant.gt_types  # extracts GTs as a numpy array
        locus = variant.CHROM + ':' + str(variant.POS)
        df_to_append = pd.DataFrame(geno, columns=[locus])  # creates a temp df to store the GTs for one locus

        # concatenate results to the main df
        df = pd.concat([df, df_to_append], axis=1)
    print("Finished subsetting VCF for window")

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

    # keep only the SNPs that are in the GWAS catalog
    gwas_snps = pd.read_csv(gwas_snp_path)
    correlation_df = correlation_df[correlation_df['locus'].isin(gwas_snps['locus'])]

    # find the SNP with the highest absolute correlation
    max_correlation_index = correlation_df['correlation'].abs().idxmax()
    max_correlation_df = correlation_df[correlation_df['locus'] == max_correlation_index]

    # add some attributes
    max_correlation_df['gene'] = gene
    max_correlation_df['str_locus'] = str_locus
    max_correlation_df['celltype'] = celltype

    # return the df as a String
    return max_correlation_df.to_csv(index=False)


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
@click.option('--coloc-dir', help='GCS file path to coloc results', type=str)
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
@click.option('--job-cpu', default=1)
@click.option('--job-storage', default='20G')
@click.command()
def main(
    snp_vcf_dir: str,
    str_vcf_dir: str,
    coloc_dir: str,
    phenotype: str,
    celltypes: str,
    gene_annotation_file: str,
    str_fdr_dir: str,
    job_cpu: int,
    job_storage: str,
):
    for celltype in celltypes.split(','):
        # read in STR eGene annotation file
        str_fdr_file = f'{str_fdr_dir}/{celltype}_qval.tsv'
        str_fdr = pd.read_csv(str_fdr_file, sep='\t')
        str_fdr = str_fdr[str_fdr['qval'] < 0.05]  # subset to eGenes passing FDR 5% threshold

        coloc_result_file = f'{coloc_dir}/{phenotype}/{celltype}/gene_summary_result.csv'
        coloc_results = pd.read_csv(coloc_result_file)
        # subset results for posterior probability of a shared causal variant >=0.5
        coloc_results = coloc_results[coloc_results['PP.H4.abf'] >= 0.5]

        # obtain inputs for LD parsing for each entry in `coloc_results`:
        for index, row in coloc_results.iterrows():
            gene = row['gene']
            # obtain snp cis-window coordinates for the gene
            gene_annotation_table = pd.read_csv(gene_annotation_file)
            gene_table = gene_annotation_table[
                gene_annotation_table['gene_ids'] == gene
            ]  # subset to particular ENSG ID
            start_snp_window = float(gene_table['start'].astype(float)) - 100000  # +-100kB window around gene
            end_snp_window = float(gene_table['end'].astype(float)) + 100000  # +-100kB window around gene
            chr = gene_table['chr'].iloc[0][3:]
            snp_window = f'{chr}:{start_snp_window}-{end_snp_window}'
            print('Obtained SNP window coordinates')

            # obtain top STR locus for the gene
            str_fdr_gene = str_fdr[str_fdr['gene_name'] == gene]
            for estr in zip(
                ast.literal_eval(str_fdr_gene['chr'].iloc[0]),
                ast.literal_eval(str_fdr_gene['pos'].iloc[0]),
            ):
                chr_num = estr[0][3:]
                pos = estr[1]
                end = str(int(pos) + 1)
                str_locus = f'{chr_num}:{pos}-{end}'
                print(f'Running LD for {gene} and {str_locus}')
                gwas_snp_path = f'{coloc_dir}/{phenotype}/{celltype}/{gene}_snp_gwas_list.csv'
                snp_vcf_path = f'{snp_vcf_dir}/chr{chr}_common_variants.vcf.bgz'
                str_vcf_path = f'{str_vcf_dir}/hail_filtered_chr{chr_num}.vcf.bgz'
                # run coloc
                b = get_batch()
                ld_job = b.new_python_job(
                    f'LD calc for {gene} and STR: {str_locus}; {celltype}',
                )
                ld_job.cpu(job_cpu)
                ld_job.storage(job_storage)
                snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_path, 'csi': snp_vcf_path + '.csi'})
                str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'csi': str_vcf_path + '.csi'})

                result = ld_job.call(
                    ld_parser,
                    snp_input,
                    str_input,
                    str_locus,
                    snp_window,
                    gwas_snp_path,
                    gene,
                    celltype,
                )

                write_path = output_path(
                    f'coloc_str_ld/{phenotype}/{celltype}/{gene}_chr{chr_num}_{pos}_ld_results.csv',
                    'analysis',
                )
                b.write_output(result.as_str(), write_path)

            b.run(wait=False)


if __name__ == '__main__':
    main()
