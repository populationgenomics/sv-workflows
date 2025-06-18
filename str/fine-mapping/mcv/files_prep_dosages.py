#!/usr/bin/env python3

"""
This script prepares the GT dosages files for input into files_prep_residualizer.py

analysis-runner --dataset tenk10k --access-level test --memory 8G --description "Dosages file prep" --output-dir str/associatr/final_freeze/fine_mapping/susie_mcv/prep_files \
files_prep_dosages.py --estrs-path=gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/cell-type-spec/estrs.csv

"""

import click
from cpg_utils.hail_batch import get_batch, output_path
import pandas as pd
from cpg_utils import to_path


def tr_extract_genotype_matrix(chrom, start, end):
    """
    Extracts the genotype matrix for tandem repeats from a VCF file.
    (summed repeats)
    """

    from cyvcf2 import VCF
    import pandas as pd

    local_file = 'local.vcf.bgz'
    local_file_index = 'local.vcf.bgz.tbi'

    tr_file = to_path(
        f'gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific/hail_filtered_{chrom}.vcf.bgz'
    )
    tr_file_index = to_path(
        f'gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific/hail_filtered_{chrom}.vcf.bgz.tbi'
    )

    tr_file.copy(local_file)
    tr_file_index.copy(local_file_index)

    vcf = VCF(local_file)

    region = f"{chrom}:{start}-{end}"
    samples = vcf.samples
    genotype_dict = {sample: {} for sample in samples}

    for variant in vcf(region):
        coord = f"{variant.CHROM}:{variant.POS}:{variant.INFO.get('RU')}"
        repcn_data = variant.format('REPCN')  # array of strings like "10/12"

        if repcn_data is None:
            continue  # skip if field missing

        for i, sample in enumerate(samples):
            val = repcn_data[i]
            try:
                if isinstance(val, bytes):
                    val = val.decode()

                if val in {"./.", ".", ""}:
                    summed = None
                else:
                    allele_1, allele_2 = map(int, val.split('/'))
                    summed = allele_1 + allele_2
            except Exception:
                summed = None

            genotype_dict[sample][coord] = summed

    df = pd.DataFrame.from_dict(genotype_dict, orient='index')
    df.index.name = "sample"
    df.sort_index(axis=1, inplace=True)
    df.reset_index(inplace=True)
    return df


def snp_extract_genotype_matrix(chrom, start, end):
    """
    Extracts the genotype matrix for SNPs from a VCF file.

    """
    from cyvcf2 import VCF
    import pandas as pd

    local_file = 'local.vcf.bgz'
    local_file_index = 'local.vcf.bgz.tbi'

    snp_file = to_path(
        f'gs://cpg-bioheart-test/tenk10k/str/associatr/common_variant_snps/hail_filtered_{chrom}.vcf.bgz'
    )
    snp_file_index = to_path(
        f'gs://cpg-bioheart-test/tenk10k/str/associatr/common_variant_snps/hail_filtered_{chrom}.vcf.bgz.tbi'
    )

    snp_file.copy(local_file)
    snp_file_index.copy(local_file_index)

    vcf = VCF(local_file)
    region = f"{chrom}:{start}-{end}"
    samples = vcf.samples
    genotype_dict = {sample: {} for sample in samples}

    for variant in vcf(region):
        coord = f"{variant.CHROM}:{variant.POS}:{variant.INFO.get('RU')}" # e.g., "chr1:123456:A-T"
        gt_data = variant.genotypes  # shape: (n_samples, 3)

        for i, sample in enumerate(samples):
            try:
                allele_1, allele_2, phased = gt_data[i]  # example: [0, 1, True]
                # Ignore phased for dosage purposes
                if allele_1 is None or allele_2 is None:
                    summed = None
                else:
                    summed = allele_1 + allele_2 #summed GT in SNP context is number of non-REF alleles
            except Exception:
                summed = None

            genotype_dict[sample][coord] = summed

    df = pd.DataFrame.from_dict(genotype_dict, orient='index')
    df.index.name = "sample"
    df.sort_index(axis=1, inplace=True)
    df.reset_index(inplace=True)
    return df


def dosages(gene):
    import pandas as pd
    import numpy as np
    from cpg_utils.hail_batch import output_path

    """
    Extracts the dosage files for each gene in the df DataFrame.
    """ ""

    # Load metadata to extract the window around the gene
    gene_info = pd.read_csv('gs://cpg-bioheart-test/tenk10k/saige-qtl/300libraries_n1925_adata_raw_var.csv')
    gene_ensg = gene
    gene_info_filtered = gene_info[gene_info['gene_ids'] == gene_ensg]
    chrom = gene_info_filtered.iloc[0]['chr']
    start = int(gene_info_filtered.iloc[0]['start'])
    end = int(gene_info_filtered.iloc[0]['end']) + 100_000
    start_body = max(0, start - 100_000)

    # Extract genotype matrices for tandem repeats and SNPs
    tr_df = tr_extract_genotype_matrix(chrom, start_body, end)
    snp_df = snp_extract_genotype_matrix(chrom, start_body, end)

    # Merge the two genotype dfs together
    variant_df = tr_df.merge(snp_df)
    variant_df['sample'] = variant_df['sample'].astype(float)

    # Save file
    variant_df.to_csv(output_path(f"dosages/{gene_ensg}_dosages.csv"), header=True)


@click.option('--estrs-path', type=str, required=True, help='Path to the estrs file.')
@click.command()
def main(estrs_path):
    b = get_batch(name='Dosages files prep for SuSie MCV')
    df = pd.read_csv(estrs_path)
    df = df.drop_duplicates(subset=['gene_name', 'chr'])
    # sort by chromosome
    for chrom in df['chr'].unique():
        df_chr = df[df['chr'] == chrom]
        for gene in df_chr['gene_name'].unique():
            if to_path(output_path(f'dosages/{gene}_dosages.csv')).exists():
                print(f"Dosage file for {gene} already exists, skipping.")
                continue
            j = b.new_python_job(name=f'Prepare for {chrom} and {gene}')
            j.cpu(1)
            j.storage('10G')
            j.call(dosages, gene)
    b.run(wait=False)


if __name__ == '__main__':
    main()
