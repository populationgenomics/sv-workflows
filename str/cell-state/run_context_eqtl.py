#!/usr/bin/env python3

"""
This script runs context eQTL analysis for a specific cell type and pathway.
It processes gene lists, extracts dosage data, and performs association tests.

analysis-runner --access-level test --dataset tenk10k --description "Run context eQTL analysis" \
--output-dir "str/cellstate/o_files/tob" run_context_eqtl.py

"""
import json
import pandas as pd
import click
from cpg_utils.hail_batch import output_path
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch

import pandas as pd

def process_gene_fast(pheno_cov_dir, gene, chromosome, cell_type, pathway):
    import pandas as pd
    import numpy as np
    import statsmodels.api as sm
    from patsy import dmatrix
    from cpg_utils.hail_batch import output_path
    from scipy.linalg import qr

    def drop_collinear_columns(X, tol=1e-10):
        X_np = np.asarray(X)  # Ensure it's a NumPy array
        Q, R, P = qr(X_np, mode='economic', pivoting=True)
        rank = np.sum(np.abs(np.diag(R)) > tol)
        keep_cols = P[:rank]
        X_reduced = X_np[:, keep_cols]
        return X_reduced, keep_cols

    # Load phenotype + covariate data
    df = pd.read_csv(f'{pheno_cov_dir}/{pathway}/{cell_type}/{chromosome}/{gene}_pheno_cov.csv')

    # Load TR dosage matrix
    dosage = pd.read_csv(f'gs://cpg-tenk10k-test/str/cellstate/input_files/dosages/{chromosome}/{gene}_dosages.csv')
    dosage['sample'] = 'CPG' + dosage['sample'].astype(str)

    # Merge by sample ID
    df = df.merge(dosage, left_on='sample_id', right_on='sample', how='inner')

    # Ensure categorical format
    df['activity'] = pd.Categorical(df['activity'], categories=['low', 'medium', 'high'], ordered=False)

    # Get outcome
    y = df['gene_inverse_normal'].values

    # Base formula (excluding genotype)
    fixed_formula = (
        "sex + age + " +
        " + ".join([f"score_{i}" for i in range(1, 13)]) + " + " +
        " + ".join([f"rna_PC{i}" for i in range(1, 7)]) + " + " +
        "C(activity) + C(sample_id)"
    )
    X_base = dmatrix(fixed_formula, df, return_type='dataframe')

    # Collect variant IDs
    variant_cols = [col for col in dosage.columns if col != "sample"]

    # Genotype matrix (N samples x M variants)
    G = df[variant_cols].copy()

    # Mean imputation for missing genotypes
    G_imputed = G.fillna(G.mean())

    # Center genotype (optional, improves interpretability & stability)
    G_centered = G_imputed - G_imputed.mean()

    # Create interaction terms with activity
    activity_dummies = pd.get_dummies(df['activity'], prefix='activity')[['activity_medium', 'activity_high']]
    interaction_med = G_centered.values * activity_dummies['activity_medium'].values[:, np.newaxis]
    interaction_high = G_centered.values * activity_dummies['activity_high'].values[:, np.newaxis]

    # Add all to design matrix
    X_all = pd.concat(
        [X_base.reset_index(drop=True),
         pd.DataFrame(G_centered.values, columns=[f'genotype_{v}' for v in variant_cols]),
         pd.DataFrame(interaction_med, columns=[f'{v}:activity_medium' for v in variant_cols]),
         pd.DataFrame(interaction_high, columns=[f'{v}:activity_high' for v in variant_cols])
         ],
        axis=1
    )

    # Drop colinear columns
    X_clean, dropped = drop_collinear_columns(X_all)

    # Fit OLS model
    model = sm.OLS(y, X_clean).fit()

    # Collect results
    results = []
    for var_id in variant_cols:
        result = {
            'gene': gene,
            'variant': var_id,
            'beta_genotype': model.params.get(f'genotype_{var_id}', np.nan),
            'pval_genotype': model.pvalues.get(f'genotype_{var_id}', np.nan),
            'beta_genox_med': model.params.get(f'{var_id}:activity_medium', np.nan),
            'pval_genox_med': model.pvalues.get(f'{var_id}:activity_medium', np.nan),
            'beta_genox_high': model.params.get(f'{var_id}:activity_high', np.nan),
            'pval_genox_high': model.pvalues.get(f'{var_id}:activity_high', np.nan),
        }
        results.append(result)

    results_df = pd.DataFrame(results)
    results_df.to_csv(
        output_path(f'context_eqtl_results/{pathway}/{cell_type}/{chromosome}/{gene}_context_eqtl_results.csv'),
        index=False,
    )

@click.option('--input-gene-list-dir', default='gs://cpg-tenk10k-test/str/cellstate/input_files/tob/scRNA_gene_lists')
@click.option('--cis-window-dir', default='gs://ccpg-tenk10k-test/str/cellstate/input_files/tob/cis_window_files')
@click.option('--pheno-cov-dir', default='gs://cpg-tenk10k-test/str/cellstate/input_files/tob/pheno_cov_csv')
@click.option('--cell-type', default='B_naive', help='Cell type to process')
@click.option('--pathway', default='GOBP_MULTI_MULTICELLULAR_ORGANISM_PROCESS_subtype', help='Pathway to process')
@click.option('--chromosome', default='chr22', help='Chromosome to process')
@click.option('--job-storage', help='Storage of the batch job eg 30G', default='4G')
@click.option('--job-memory', help='Memory of the batch job', default='standard')
@click.option('--job-cpu', help='Number of CPUs of Hail batch job', default=0.25)
@click.command()
def main(
    input_gene_list_dir, cis_window_dir, pheno_cov_dir, cell_type, pathway, chromosome, job_storage, job_memory, job_cpu
):
    """
    Run context eQTL analysis for a specific cell type and pathway
    """
    b = get_batch(name='run context eQTL analysis')
    gene_list_path = (
        input_gene_list_dir
        + f'/{pathway}/1_min_pct_cells_expressed/{cell_type}/{chromosome}_{cell_type}_gene_list.json'
    )
    with open(to_path(gene_list_path), 'r') as f:
        gene_list = json.load(f)

    # Process each gene in the list
    for gene in gene_list:
        if to_path(
            output_path(f'context_eqtl_results/{pathway}/{cell_type}/{chromosome}/{gene}_context_eqtl_results.csv')
        ).exists():
            print(f"Context eQTL results for {gene} already exist, skipping.")
            #continue
        j = b.new_python_job(
            name=f'{cell_type}: {chromosome} {gene}',
        )
        j.cpu(job_cpu)
        j.memory(job_memory)
        j.storage(job_storage)
        j.call(process_gene_fast, pheno_cov_dir, gene, chromosome, cell_type, pathway)

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter