#!/usr/bin/env python3

"""
This script prepares the X and Y residualized files for SuSiE multiple causal variant assumption mapping.

analysis-runner --dataset bioheart --access-level test --description "Residualized files prep for SuSie MCV" --output-dir tenk10k/str/associatr/final_freeze/fine_mapping/susie_mcv/prep_files \
files_prep.py --estrs-path=gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/bioheart_n975_and_tob_n950/TableS1.csv
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


def process_cohort(variant_df, ycov, gene_ensg):
    """
    Aligns phenotype/covariate numpy array with variant_df,
    residualizes both y and X with respect to covariates,
    and returns residualized phenotype and genotype matrix.

    Parameters:
    - variant_df: pandas DataFrame with 'sample' column and variant columns
    - ycov: numpy array with columns: sample_id, phenotype, covariates
    - gene_ensg: string (used only for naming, can be optional)

    Returns:
    - y_resid_series: pandas Series (sample-indexed residualized phenotype)
    - X_resid_df: pandas DataFrame (sample-indexed residualized genotypes)
    """
    import pandas as pd
    import numpy as np
    import statsmodels.api as sm

    # Define covariate column names
    n_covariates = ycov.shape[1] - 2
    ycov_columns = ['sample', 'phenotype'] + [f'covar{i+1}' for i in range(n_covariates)]

    # Convert to DataFrame
    ycov_df = pd.DataFrame(ycov, columns=ycov_columns)
    ycov_df['phenotype'] = ycov_df['phenotype'].astype(float)
    for c in ycov_columns[2:]:
        ycov_df[c] = ycov_df[c].astype(float)

    # Merge on 'sample'
    merged = pd.merge(variant_df, ycov_df, on='sample', how='inner')

    # Extract genotype matrix and covariates
    variant_cols = variant_df.columns.drop('sample')
    X = merged[variant_cols].copy()
    X_imputed = X.apply(lambda col: col.fillna(col.mean()), axis=0).values

    y = merged['phenotype'].values
    C = merged[ycov_columns[2:]].values  # only covariates
    C = sm.add_constant(C)

    # Residualize y
    model_y = sm.OLS(y, C).fit()
    y_resid = model_y.resid

    # Residualize X
    X_resid = np.empty_like(X_imputed)
    for j in range(X_imputed.shape[1]):
        model = sm.OLS(X_imputed[:, j], C).fit()
        X_resid[:, j] = model.resid

    # Format outputs
    y_resid_series = pd.Series(y_resid, index=merged['sample'], name='phenotype_resid')
    X_resid_df = pd.DataFrame(X_resid, columns=variant_cols, index=merged['sample'])

    return y_resid_series, X_resid_df


def residualizer(cell_type, gene_name, snp_input, tr_input):
    from cpg_utils import to_path
    import pandas as pd
    import numpy as np
    from cpg_utils.hail_batch import output_path

    # === LOAD metadata ===
    gene_info = pd.read_csv('gs://cpg-bioheart-test/tenk10k/saige-qtl/300libraries_n1925_adata_raw_var.csv')
    gene_ensg = gene_name
    gene_info_filtered = gene_info[gene_info['gene_ids'] == gene_ensg]
    chrom = gene_info_filtered.iloc[0]['chr']
    start = int(gene_info_filtered.iloc[0]['start'])
    end = int(gene_info_filtered.iloc[0]['end']) + 100_000
    start_body = max(0, start - 100_000)

    tr_df = tr_extract_genotype_matrix(tr_input['vcf'], chrom, start_body, end)
    snp_df = snp_extract_genotype_matrix(snp_input['vcf'], chrom, start_body, end)
    variant_df = tr_df.merge(snp_df)
    variant_df['sample'] = variant_df['sample'].astype(float)

    # === REMOVE INDELS THAT LOOK LIKE TRs === #
    meta = pd.read_csv(
        f'gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/snps_and_strs/bioheart_n975_and_tob_n950/rm_str_indels_dup_strs/meta_results/{cell_type}/{chrom}/{gene_ensg}_100000bp_meta_results.tsv',
        sep='\t',
    )
    meta['variant_id'] = meta['chr'].astype(str) + ':' + meta['pos'].astype(str) + ':' + meta['motif'].astype(str)
    columns_to_keep = ['sample'] + [col for col in variant_df.columns if col in meta['variant_id'].values]
    variant_df = variant_df[columns_to_keep]

    # === COHORT 1 ===
    ycov_1 = np.load(
        to_path(
            f'gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/bioheart_n975/pheno_cov_numpy/v1/{cell_type}/{chrom}/{gene_ensg}_pheno_cov.npy'
        )
    )
    y_resid_1, X_resid_1 = process_cohort(variant_df, ycov_1, gene_ensg)

    # === COHORT 2 ===
    ycov_2 = np.load(
        to_path(
            f'gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/tob_n950/pheno_cov_numpy/v1/{cell_type}/{chrom}/{gene_ensg}_pheno_cov.npy'
        )
    )
    y_resid_2, X_resid_2 = process_cohort(variant_df, ycov_2, gene_ensg)

    # === STACK & EXPORT ===
    y_resid_combined = pd.concat([y_resid_1, y_resid_2])
    X_resid_combined = pd.concat([X_resid_1, X_resid_2])

    # Sort by sample ID just in case
    y_resid_combined = y_resid_combined.sort_index()
    X_resid_combined = X_resid_combined.sort_index()

    # Save files
    y_resid_combined.to_csv(output_path(f"{gene_ensg}_{cell_type}_meta_cleaned_y_resid.csv"), header=True)
    X_resid_combined.to_csv(output_path(f"{gene_ensg}_{cell_type}_meta_cleaned_X_resid.csv"))


@click.option('--estrs-path', type=str, required=True, help='Path to the estrs file.')
@click.command()
def main(estrs_path):
    b = get_batch(name='Residualized files prep for SuSie MCV')
    df = pd.read_csv(estrs_path, sep='\t', index_col=0)
    df = df.drop_duplicates(subset=['cell_type', 'gene_name'])
    # sort by chromosome
    for chrom in df['chr'].unique():
        df_chr = df[df['chr'] == chrom]
        tr_vcf_file = to_path(
            f'gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific/hail_filtered_{chrom}.vcf.bgz'
        )
        snp_vcf_file = to_path(
            f'gs://cpg-bioheart-test/tenk10k/str/associatr/common_variant_snps/hail_filtered_{chrom}.vcf.bgz'
        )
        snp_input = get_batch().read_input_group(**{'vcf': snp_vcf_file, 'tbi': snp_vcf_file + '.tbi'})
        tr_input = get_batch().read_input_group(**{'vcf': tr_vcf_file, 'tbi': tr_vcf_file + '.tbi'})
        for row in df_chr.itertuples(index=False):
            j = b.new_python_job(f'Prepare for {row.cell_type} {row.gene_name}')
            j.cpu(0.25)
            j.storage('5G')
            j.call(residualizer, row.cell_type, row.gene_name, snp_input, tr_input)
    b.runs(wait=False)


if __name__ == '__main__':
    main()
