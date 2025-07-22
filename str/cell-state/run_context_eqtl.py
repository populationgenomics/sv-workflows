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
    import os
    from cpg_utils.hail_batch import output_path

    # Load phenotype + covariate data
    df = pd.read_csv(f'{pheno_cov_dir}/{pathway}/{cell_type}/{chromosome}/{gene}_pheno_cov.csv')

    # Load genotype matrix
    dosage = pd.read_csv(f'gs://cpg-tenk10k-test/str/cellstate/input_files/dosages/{chromosome}/{gene}_dosages.csv')
    dosage['sample'] = 'CPG' + dosage['sample'].astype(str)

    # Merge
    df = df.merge(dosage, left_on='sample_id', right_on='sample', how='inner')

    # Clean categories
    df['activity'] = pd.Categorical(df['activity'], categories=['low', 'medium', 'high'], ordered=False)

    # Prepare fixed covariates
    covariate_cols = ['sex', 'age'] + [f'score_{i}' for i in range(1, 13)] + [f'rna_PC{i}' for i in range(1, 7)]
    activity_dummies = pd.get_dummies(df['activity'], prefix='activity')[['activity_medium', 'activity_high']]
    sample_dummies = pd.get_dummies(df['sample_id'], prefix='sample')

    # Combine fixed covariates
    X_fixed = pd.concat([df[covariate_cols], activity_dummies, sample_dummies], axis=1).astype(float).values
    y = df['gene_inverse_normal'].values

    # Genotype variants
    variant_cols = [col for col in dosage.columns if col != "sample"]
    G = df[variant_cols].copy().fillna(df[variant_cols].mean()).astype(float)
    G -= G.mean()  # Center

    # Fit OLS per variant
    results = []
    for v in variant_cols:
        g = G[v].values
        gx_med = g * activity_dummies['activity_medium'].values
        gx_high = g * activity_dummies['activity_high'].values

        # Combine into design matrix
        X = np.column_stack([X_fixed, g, gx_med, gx_high])

        # OLS via least squares
        beta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
        df_resid = len(y) - rank
        y_hat = X @ beta
        rss = np.sum((y - y_hat) ** 2)
        mse = rss / df_resid
        se = np.sqrt(mse * np.diag(np.linalg.pinv(X.T @ X)))

        # Extract coefficients and p-values
        beta_g, beta_gx_med, beta_gx_high = beta[-3:]
        se_g, se_gx_med, se_gx_high = se[-3:]
        from scipy.stats import t
        tvals = np.array([beta_g, beta_gx_med, beta_gx_high]) / np.array([se_g, se_gx_med, se_gx_high])
        pvals = 2 * t.sf(np.abs(tvals), df_resid)

        results.append({
            'gene': gene,
            'variant': v,
            'beta_genotype': beta_g,
            'pval_genotype': pvals[0],
            'beta_genox_med': beta_gx_med,
            'pval_genox_med': pvals[1],
            'beta_genox_high': beta_gx_high,
            'pval_genox_high': pvals[2],
        })

    # Save
    pd.DataFrame(results).to_csv(
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