#!/usr/bin/env python3

"""
This script prepares the X and Y residualized files for SuSiE multiple causal variant assumption mapping.

analysis-runner --dataset tenk10k --access-level test --memory 8G --description "Residualized files prep for SuSie MCV" --output-dir str/associatr/final_freeze/fine_mapping/susie_mcv/prep_files \
files_prep_residualizer.py --estrs-path=gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/cell-type-spec/estrs.csv
"""

import click
from cpg_utils.hail_batch import get_batch, output_path
import pandas as pd
from cpg_utils import to_path


def process_cohort_wide(variant_df, ycov1, ycov2, cell_type):
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

    # STEP 1: Create a df based on pseudobulk phenotype values; and joint cohort genotype PCs and RNA PCs; age and sex

    # pull out sample_id and phenotype (pseudobulk value from ycov) from TOB and BioHEART each (dimensions are the same)

    # start with the first cohort.
    n_covariates = ycov1.shape[1] - 2
    ycov_columns = ['sample_id', 'phenotype'] + [f'covar{i+1}' for i in range(n_covariates)]
    ycov_df_1 = pd.DataFrame(ycov1, columns=ycov_columns)
    ycov_df_1 = ycov_df_1[['sample_id', 'phenotype']]
    # repeat for the second cohort
    ycov_df_2 = pd.DataFrame(ycov2, columns=ycov_columns)
    ycov_df_2 = ycov_df_2[['sample_id', 'phenotype']]
    # concatenate the two cohorts
    ycov_df = pd.concat([ycov_df_1, ycov_df_2], ignore_index=True)
    ycov_df['sample_id'] = 'CPG' + ycov_df['sample_id'].astype(int).astype(str)

    # pull out sample_id, age, sex, and first 12 geno PCs
    wide_pcs = pd.read_csv(
        'gs://cpg-tenk10k-test-analysis/saige-qtl/tenk10k-genome-2-3-eur/input_files/241210/covariates/sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv'
    )
    wide_pcs = wide_pcs.iloc[:, :15]

    # pull out sample_id and first 6 RNA PCs of the relevant clel type
    rna_pcs = pd.read_csv(
        f'gs://cpg-tenk10k-test/str/pseudobulk_finemap_mcv/n1925/rna_pcs/covariates/10_rna_pcs/{cell_type}_covariates.csv'
    )
    rna_pcs = rna_pcs.iloc[:, :7]

    # sequential merging
    ycov_df = ycov_df.merge(wide_pcs)
    ycov_df = ycov_df.merge(rna_pcs)
    ycov_df = ycov_df.rename(columns={'sample_id': 'sample'})

    # --------------------
    # STEP 2: Merge on sample_id to align
    # --------------------
    merged = pd.merge(variant_df, ycov_df, on='sample', how='inner')

    # --------------------
    # STEP 3: Extract and mean-impute X (missing GTs will not work.)
    # --------------------
    variant_cols = variant_df.columns.drop('sample')
    X = merged[variant_cols].copy()
    X_imputed = X.apply(lambda col: col.fillna(col.mean()), axis=0).values

    # Extract y and covariates
    y = merged['phenotype'].values
    cols = [col for col in ycov_df.columns if col not in ['sample', 'phenotype', 'sample_id']]

    C = merged[cols].values  # covariates only

    # --------------------
    # STEP 4: Residualize y and X
    # --------------------
    C = sm.add_constant(C)

    # Residualize phenotype
    model_y = sm.OLS(y, C).fit()
    y_resid = model_y.resid

    # Residualize each variant (the X's)
    X_resid = np.empty_like(X_imputed)
    for j in range(X_imputed.shape[1]):
        model = sm.OLS(X_imputed[:, j], C).fit()
        X_resid[:, j] = model.resid

    # --------------------
    # STEP 5: Output residualized data
    # --------------------
    y_resid_series = pd.Series(y_resid, index=merged['sample'], name='phenotype_resid')
    X_resid_df = pd.DataFrame(X_resid, columns=variant_cols, index=merged['sample'])

    return y_resid_series, X_resid_df


def residualizer(gene_name, cell_type):
    from cpg_utils import to_path
    import pandas as pd
    import numpy as np
    from cpg_utils.hail_batch import output_path

    # === LOAD metadata ===
    gene_info = pd.read_csv('gs://cpg-bioheart-test/tenk10k/saige-qtl/300libraries_n1925_adata_raw_var.csv')

    gene_ensg = gene_name
    gene_info_filtered = gene_info[gene_info['gene_ids'] == gene_ensg]
    chrom = gene_info_filtered.iloc[0]['chr']

    # === LOAD VARIANT DOSAGES ===
    variant_df = pd.read_csv(
        f'gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/fine_mapping/susie_mcv/prep_files/dosages/{gene_ensg}_dosages.csv'
    )
    variant_df['sample'] = 'CPG' + variant_df['sample'].astype(int).astype(str)

    # === REMOVE INDELS THAT LOOK LIKE TRs === # (keeping TR-looking indels can mess up fine-mapping when comparing their signal with the TR)
    meta = pd.read_csv(
        f'gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/trs_snps/rm_str_indels_dup_strs_v3/{cell_type}/{chrom}/{gene_ensg}_100000bp_meta_results.tsv',
        sep='\t',
    )  # a bit of a hack - this file contains summary stats for a list of variants where TR-looking indels and duplicate TRs have been removed.
    meta['variant_id'] = meta['chr'].astype(str) + ':' + meta['pos'].astype(str) + ':' + meta['motif'].astype(str)
    columns_to_keep = ['sample'] + [col for col in variant_df.columns if col in meta['variant_id'].values]
    variant_df = variant_df[columns_to_keep]

    # === Load in bioheart data ===
    ycov_1 = np.load(
        to_path(
            f'gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/bioheart_n975/pheno_cov_numpy/v1/{cell_type}/{chrom}/{gene_ensg}_pheno_cov.npy'
        )
    )
    # === Load in TOB data ===
    ycov_2 = np.load(
        to_path(
            f'gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/tob_n950/pheno_cov_numpy/v1/{cell_type}/{chrom}/{gene_ensg}_pheno_cov.npy'
        )
    )
    y_resid_combined, X_resid_combined = process_cohort_wide(variant_df, ycov_1, ycov_2, cell_type)

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
    df = pd.read_csv(estrs_path)
    df = df.drop_duplicates(subset=['cell_type', 'gene_name', 'chr'])
    # sort by chromosome
    for chrom in df['chr'].unique():
        df_chr = df[df['chr'] == chrom]
        for cell_type in df_chr['cell_type'].unique():
            df_cell = df_chr[df_chr['cell_type'] == cell_type]
            for gene in df_cell['gene_name'].unique():
                if to_path(output_path(f"{gene}_{cell_type}_meta_cleaned_y_resid.csv")).exists():
                    continue
                j = b.new_python_job(f'Prepare for {cell_type} {chrom}: {gene}')
                j.cpu(0.25)
                j.storage('5G')
                j.call(residualizer, gene, cell_type)
    b.run(wait=False)


if __name__ == '__main__':
    main()
