#!/usr/bin/env python3

"""
This script prepares the GT dosages files for input into files_prep_residualizer.py

analysis-runner --dataset bioheart --access-level test --memory 8G --description "Residualized files prep for SuSie MCV" --output-dir tenk10k/str/associatr/final_freeze/fine_mapping/susie_mcv/prep_files \
files_prep_dosages.py --estrs-path=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/TableS1.csv

"""

import click
from cpg_utils.hail_batch import get_batch, output_path
import pandas as pd
from cpg_utils import to_path


def tr_extract_genotype_matrix(vcf_path, chrom, start, end):
    """
    Extracts the genotype matrix for tandem repeats from a VCF file.
    (summed repeats)
    """

    from cyvcf2 import VCF
    import pandas as pd

    vcf = VCF(vcf_path)
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


def snp_extract_genotype_matrix(vcf_path, chrom, start, end):
    """
    Extracts the genotype matrix for SNPs from a VCF file.

    """
    from cyvcf2 import VCF
    import pandas as pd

    vcf = VCF(vcf_path)
    region = f"{chrom}:{start}-{end}"
    samples = vcf.samples
    genotype_dict = {sample: {} for sample in samples}

    for variant in vcf(region):
        coord = f"{variant.CHROM}:{variant.POS}:{variant.INFO.get('RU') or variant.REF + '/' + ','.join(variant.ALT)}"
        gt_data = variant.genotypes  # shape: (n_samples, 3)

        for i, sample in enumerate(samples):
            try:
                allele_1, allele_2, phased = gt_data[i]  # example: [0, 1, True]
                # Ignore phased for dosage purposes
                if allele_1 is None or allele_2 is None:
                    summed = None
                else:
                    summed = allele_1 + allele_2
            except Exception:
                summed = None

            genotype_dict[sample][coord] = summed

    df = pd.DataFrame.from_dict(genotype_dict, orient='index')
    df.index.name = "sample"
    df.sort_index(axis=1, inplace=True)
    df.reset_index(inplace=True)
    return df


def dosages(gene, snp_input, tr_input):
    import pandas as pd
    import numpy as np
    from cpg_utils.hail_batch import output_path

    """
    Extracts the dosage files for each gene in the df DataFrame.
    """ ""

    # === LOAD metadata ===
    gene_info = pd.read_csv('gs://cpg-bioheart-test/tenk10k/saige-qtl/300libraries_n1925_adata_raw_var.csv')

    gene_ensg = gene
    gene_info_filtered = gene_info[gene_info['gene_ids'] == gene_ensg]
    chrom = gene_info_filtered.iloc[0]['chr']
    start = int(gene_info_filtered.iloc[0]['start'])
    end = int(gene_info_filtered.iloc[0]['end']) + 100_000
    start_body = max(0, start - 100_000)

    tr_df = tr_extract_genotype_matrix(tr_input['vcf'], chrom, start_body, end)
    snp_df = snp_extract_genotype_matrix(snp_input['vcf'], chrom, start_body, end)
    variant_df = tr_df.merge(snp_df)
    variant_df['sample'] = variant_df['sample'].astype(float)

    # Save files
    variant_df.to_csv(output_path(f"dosages/{gene_ensg}_dosages.csv"), header=True)


@click.option('--estrs-path', type=str, required=True, help='Path to the estrs file.')
@click.command()
def main(estrs_path):
    b = get_batch(name='Residualized files prep for SuSie MCV')
    df = pd.read_csv(estrs_path)
    df = df.drop_duplicates(subset=['gene_name', 'chr'])
    # sort by chromosome
    for chrom in df['chr'].unique():
        df_chr = df[df['chr'] == chrom]
        for gene in df_chr['gene_name'].unique():
            tr_vcf_file = f'gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific/hail_filtered_{chrom}.vcf.bgz'
            snp_vcf_file = (
                f'gs://cpg-bioheart-test/tenk10k/str/associatr/common_variant_snps/hail_filtered_{chrom}.vcf.bgz'
            )

            snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_file, 'tbi': snp_vcf_file + '.tbi'})
            tr_input = get_batch().read_input_group(**{'vcf': tr_vcf_file, 'tbi': tr_vcf_file + '.tbi'})

            j = b.new_python_job(name=f'Prepare for {chr} and {gene}')
            j.cpu(0.25)
            j.storage('5G')
            j.call(dosages, gene, snp_input, tr_input)
        break  # try only one chromosome for now
    b.run(wait=False)


if __name__ == '__main__':
    main()
