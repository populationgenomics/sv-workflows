#!/usr/bin/env python3

"""

This script calculates pairwise LD (Pearson correlation) between fine-mapped eTRs and variants in the GWAS catalog.

Workflow:
1) Subset GWAS catalog to the cis window.
2) Extract the genotypes for the lead eSTR and other variants in the cis window (and in the GWAS catalog).
3) Calculate pairwise correlation of the lead eSTR locus and other variants.
4) Save results of the top correlated variant in a TSV file.

analysis-runner --dataset "bioheart" \
    --description "Calculate LD for fine-mapped eSTRs with GWAS variants" \
    --access-level "test" \
    --memory='8G' \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr/coloc-ld/fm_strs_only/v4" \
    coloc_ld_fm_parser.py \
    --fm-csv=gs://cpg-bioheart-test/str/associatr/coloc/estrs_fm_coloc_list_for_ld.csv


"""


import click
import pandas as pd

from cpg_utils.hail_batch import get_batch
from cpg_utils import to_path


def ld_parser(
    pheno_csv,
    fm_chrom,
    snp_input,
    str_input,
):
    import numpy as np
    from cyvcf2 import VCF

    max_corr_master_df = pd.DataFrame()
    pheno_df = pd.read_csv(pheno_csv, sep='\t')
    gene_annotation_table = pd.read_csv(
        'gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv',
    )
    for index,row in fm_chrom.iterrows():
        chrom = row['chr']
        pos = row['pos']
        end = row['end']
        motif = row['motif_x']
        pheno = row['pheno']
        gene = row['gene_name']
        lead_str_locus = f'{chrom}_{pos}_{motif}'

        # Extract the cis window for the gene
        gene_table = gene_annotation_table[gene_annotation_table['gene_ids'] == gene]  # subset to particular ENSG ID
        start_window = float(gene_table['start'].astype(float)) - 100000  # +-100kB window around gene
        end_window = float(gene_table['end'].astype(float)) + 100000  # +-100kB window around gene
        chr = gene_table['chr'].iloc[0]
        cis_window = f'{chr}:{start_window}-{end_window}'
        print('Obtained cis_window coordinates')

        # Subset the GWAS catalog to the cis window
        pheno_df_cis = pheno_df[
            (pheno_df['chromosome'] == chrom)
            & (pheno_df['position'] >= start_window)
            & (pheno_df['position'] <= end_window)
        ]
        print('Subsetted GWAS catalog to the cis window')

        # Extract the genotypes for SNVs in the cis window
        snp_df = pd.DataFrame(columns=['individual'])
        snp_vcf = VCF(snp_input['vcf'])
        snp_df['individual'] = snp_vcf.samples

        print('Starting to subset SNP VCF for window...')
        for variant in snp_vcf(cis_window):
            gt = variant.gt_types  # extracts GTs as a numpy array
            gt[gt == 3] = 2
            motif = str(variant.INFO.get('RU')).replace('-', '_')
            snp = variant.CHROM + '_' + str(variant.POS) + '_' + motif
            df_to_append = pd.DataFrame(gt, columns=[snp])  # creates a temp df to store the GTs for one locus
            snp_df = pd.concat([snp_df, df_to_append], axis=1)
        print('Finished reading SNP VCF')

        # Extract the genotypes for STRs in the cis window
        str_df = pd.DataFrame(columns=['individual'])
        str_vcf = VCF(str_input['vcf'])
        str_df['individual'] = str_vcf.samples
        print('Starting to subset STR VCF for window...')
        for variant in str_vcf(cis_window):
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
            snp = variant.CHROM + '_' + str(variant.POS) + '_' + str(variant.INFO.get('RU'))

            str_df[snp] = sums

        # merge the STR and SNP GT dfs together
        merged_df = pd.merge(str_df, snp_df, on='individual')

        # calculate pairwise LD with the lead STR
        correlation_series = merged_df.drop(columns='individual').corrwith(merged_df[lead_str_locus])
        correlation_df = pd.DataFrame(correlation_series, columns=['correlation'])
        correlation_df['locus'] = correlation_df.index
        # keep only the variants in the GWAS catalog
        correlation_df = correlation_df[correlation_df['locus'].isin(pheno_df_cis['snp'])]

        # find the VARIANT with the highest absolute correlation
        max_correlation_index = correlation_df['correlation'].abs().idxmax()
        max_correlation_row = correlation_df[correlation_df['locus'] == max_correlation_index]

        # add some attributes
        max_correlation_row['gene'] = gene
        max_correlation_row['lead_str_locus'] = lead_str_locus
        max_correlation_row['max_corr_gwas_locus'] = max_correlation_index
        max_correlation_row['pheno'] = pheno
        max_correlation_row = max_correlation_row[['gene', 'lead_str_locus', 'pheno', 'correlation', 'max_corr_gwas_locus']]

        max_corr_master_df = pd.concat([max_corr_master_df, max_correlation_row], axis=0)

    max_corr_master_df.to_csv(
        f'gs://cpg-bioheart-test-analysis/str/associatr/coloc-ld/fm_strs_only/v4/{pheno}/{chrom}/{pheno}_{chrom}_corr.tsv',
        sep='\t',
        index=False,
    )


@click.option('--snp-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/tob_freeze_1/bgzip_tabix/v4')
@click.option('--str-vcf-dir', default='gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1-chr-specific')
@click.option('--fm-csv', required=True, help='Fine-mapped eSTRs CSV file path')
@click.command()
def main(fm_csv, snp_vcf_dir, str_vcf_dir):
    b = get_batch(name='Calculate LD for fine-mapped eSTRs with GWAS variants')
    fm = pd.read_csv(fm_csv)
    fm = fm.drop_duplicates(subset=['chr', 'pos', 'end', 'motif_x', 'pheno'])
    pheno_list = fm['pheno'].unique()
    # map pheno to the pheno csv file path:

    for pheno in pheno_list:
        fm_pheno = fm[fm['pheno'] == pheno]
        for chrom in fm_pheno['chr'].unique():
            if to_path(f'gs://cpg-bioheart-test-analysis/str/associatr/coloc-ld/fm_strs_only/v4/{pheno}/{chrom}/{pheno}_{chrom}_corr.tsv').exists():
                print(f'File already exists for {pheno} and {chrom}')
                continue
            pheno_map = {
            'sle_GCST003156': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/bentham_2015_26502338_sle_parsed.tsv',
            'ra_GCST90132223': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90132223_parsed.tsv',
            'covid_GCST011071': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv',
            'ibd_liu2023': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/ibd_EAS_EUR_SiKJEF_meta_IBD.tsv',
            'alzheimer_GCST90027158': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90027158.h_parsed.tsv',
            'breastca_GCST004988': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004988.h_parsed.tsv',
            'colorectalca_GCST90129505': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/colorectalca_GCST90129505_parsed.tsv',
            'lungca_GCST004748': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST004748.h_parsed.tsv',
            'lymphoma_GCST90018878': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90018878.h_parsed.tsv',
            'parkinson_GCST009325': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST009325.h_parsed.tsv',
            'prostateca_GCST90274713': 'gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST90274713.h_parsed.tsv',
            'gymrek-ukbb-alanine-aminotransferase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_alanine_aminotransferase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-albumin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_albumin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-alkaline_phosphatase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_alkaline_phosphatase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-apolipoprotein_a': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_apolipoprotein_a_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-apolipoprotein_b': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_apolipoprotein_b_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-aspartate_aminotransferase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_aspartate_aminotransferase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-c_reactive_protein': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_c_reactive_protein_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-calcium': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_calcium_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-cholesterol': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_cholesterol_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-creatinine': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_creatinine_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-cystatin_c': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_cystatin_c_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-eosinophil_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_eosinophil_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-eosinophil_percent': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_eosinophil_percent_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-gamma_glutamyltransferase': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_gamma_glutamyltransferase_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-glucose': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_glucose_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-glycated_haemoglobin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_glycated_haemoglobin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-haematocrit': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_haematocrit_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-haemoglobin_concentration': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_haemoglobin_concentration_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-hdl_cholesterol': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_hdl_cholesterol_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-igf_1': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_igf_1_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-ldl_cholesterol_direct': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_ldl_cholesterol_direct_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-phosphate': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_phosphate_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-shbg': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_shbg_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-total_bilirubin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_total_bilirubin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-total_protein': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_total_protein_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-triglycerides': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_triglycerides_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-urate': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_urate_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-urea': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_urea_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-vitamin_d': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_vitamin_d_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-white_blood_cell_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_white_blood_cell_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-lymphocyte_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_lymphocyte_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-lymphocyte_percent': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_lymphocyte_percent_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_corpuscular_haemoglobin': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_corpuscular_haemoglobin_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_corpuscular_haemoglobin_concentration': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_corpuscular_haemoglobin_concentration_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_corpuscular_volume': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_corpuscular_volume_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_platelet_volume': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_platelet_volume_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-mean_sphered_cell_volume': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_mean_sphered_cell_volume_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-neutrophil_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_neutrophil_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-neutrophil_percent': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_neutrophil_percent_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-platelet_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_platelet_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-platelet_crit': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_platelet_crit_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-platelet_distribution_width': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_platelet_distribution_width_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-red_blood_cell_count': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_red_blood_cell_count_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            'gymrek-ukbb-red_blood_cell_distribution_width': f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_red_blood_cell_distribution_width_snp_str_gwas_results_hg38_{chrom}.tab.gz',
            }
            ld_job = b.new_python_job(
                f'LD calc for {chrom}; {pheno}',
            )
            ld_job.cpu(1)
            ld_job.storage('10G')


            fm_pheno_chrom = fm_pheno[fm_pheno['chr'] == chrom]
            pheno_csv = pheno_map[pheno]

            snp_vcf_path = f'{snp_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'
            str_vcf_path = f'{str_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'

            snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_path, 'tbi': snp_vcf_path + '.tbi'})
            str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'})

            ld_job.call(
                ld_parser,
                pheno_csv,
                fm_pheno_chrom,
                snp_input,
                str_input,
            )

    b.run(wait=False)


if __name__ == '__main__':
    main()
