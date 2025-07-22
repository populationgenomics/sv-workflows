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

    """
    Efficient matrix-based regression of genotype × activity interaction using mean-imputed genotype data.

    Args:
        pheno_cov_dir (str): Directory with phenotype + covariate data.
        gene (str): Ensembl or gene symbol.
        chromosome (str): Chromosome name (e.g., 'chr22').
        cell_type (str): Cell type name.
        pathway (str): Pathway identifier (used as subfolder in GCS).
    """

    # Load phenotype + covariate data
    df = pd.read_csv(f'{pheno_cov_dir}/{pathway}/{cell_type}/{chromosome}/{gene}_pheno_cov.csv')

    # Load TR dosage matrix
    dosage = pd.read_csv(f'gs://cpg-tenk10k-test/str/cellstate/input_files/dosages/{chromosome}/{gene}_dosages.csv')
    dosage['sample'] = 'CPG' + dosage['sample'].astype(str)

    # Merge by sample ID
    df = df.merge(dosage, left_on='sample_id', right_on='sample', how='inner')

    # Ensure categorical encoding of activity_bin
    df['activity'] = pd.Categorical(df['activity'], categories=['low', 'medium', 'high'], ordered=False)

    # Outcome vector
    y = df['gene_inverse_normal'].values

    # Precompute fixed covariate matrix (excluding genotype terms)
    fixed_formula = (
        "sex + age + " +
        " + ".join([f"score_{i}" for i in range(1, 13)]) + " + " +
        " + ".join([f"rna_PC{i}" for i in range(1, 7)]) + " + " +
        "C(activity) + C(sample_id)"
    )
    X_fixed = dmatrix(fixed_formula, df, return_type='dataframe')

    # Identify genotype variant columns
    variant_cols = [col for col in dosage.columns if col != "sample"]

    # Extract and mean-impute genotype matrix
    G_raw = df[variant_cols].copy()
    G_imputed = G_raw.fillna(G_raw.mean())  # mean imputation
    G_centered = G_imputed - G_imputed.mean()  # center genotypes

    # Build interaction terms
    interaction_terms = {}
    for level in ['medium', 'high']:
        interaction_terms[f'genotype:activity_{level}'] = (
            G_centered.values * (df['activity'] == level).values[:, None]
        )

    # Assemble design matrix
    X_full = pd.concat(
        [X_fixed.reset_index(drop=True),
         pd.DataFrame(G_centered.values, columns=[f'G_{v}' for v in variant_cols]),
         pd.DataFrame(interaction_terms['genotype:activity_medium'], columns=[f'G_{v}:activity_medium' for v in variant_cols]),
         pd.DataFrame(interaction_terms['genotype:activity_high'], columns=[f'G_{v}:activity_high' for v in variant_cols])],
        axis=1
    )

    # Fit model for all variants jointly using matrix regression
    model = sm.OLS(y, X_full).fit()

    # Extract coefficients and p-values for genotype and interaction terms
    results = []
    for var in variant_cols:
        results.append({
            "gene": gene,
            "variant": var,
            "beta_genotype": model.params.get(f'G_{var}', np.nan),
            "pval_genotype": model.pvalues.get(f'G_{var}', np.nan),
            "beta_genox_med": model.params.get(f'G_{var}:activity_medium', np.nan),
            "pval_genox_med": model.pvalues.get(f'G_{var}:activity_medium', np.nan),
            "beta_genox_high": model.params.get(f'G_{var}:activity_high', np.nan),
            "pval_genox_high": model.pvalues.get(f'G_{var}:activity_high', np.nan),
        })

    # Save results
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