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

def process_gene(pheno_cov_dir, gene, chromosome, cell_type, pathway):
    from cpg_utils.hail_batch import output_path
    import statsmodels.formula.api as smf
    import pandas as pd
    import numpy as np

    """
    Run genotype × activity_bin interaction model for each TR variant near a gene.

    Args:
        pseudobulk_dir (str): Directory containing per-gene pseudobulk phenotype + covariate data.
        cis_window_dir (str): Not used here, but kept for possible BED file integration.
        gene (str): Gene symbol.
        chromosome (str): Chromosome (e.g., "chr22").
        cell_type (str): Cell type name.
        pathway (str): Pathway identifier.

    Returns:
        results (list of dict): Association statistics for each tested variant.
    """
    # Load phenotype + covariate data
    df = pd.read_csv(f'{pheno_cov_dir}/{pathway}/{cell_type}/{chromosome}/{gene}_pheno_cov.csv')

    # Load dosage data for all TRs near this gene
    dosage = pd.read_csv(f'gs://cpg-tenk10k-test/str/cellstate/input_files/dosages/{chromosome}/{gene}_dosages.csv')
    dosage['sample'] = 'CPG' + dosage['sample'].astype(str)

    # Merge dosage into phenotype data
    df = df.merge(dosage, left_on='sample_id', right_on='sample', how='inner')

    # Make sure activity_bin is categorical
    df["activity"] = pd.Categorical(df["activity"], categories=["low", "medium", "high"], ordered=False)

    results = []

    # Iterate over TR variants
    variant_cols = [col for col in dosage.columns if col != "sample"]

    for var_id in variant_cols:
        df["genotype"] = df[var_id]

        try:
            df_model = df.dropna(subset=["genotype", "gene_inverse_normal"])
            model = smf.ols(
                "gene_inverse_normal ~ genotype * C(activity) + sex+age+score_1+score_2+score_3+score_4+score_5+score_6+score_7+score_8+score_9+score_10+score_11+score_12+rna_PC1+rna_PC2+rna_PC3+rna_PC4+rna_PC5+rna_PC6+ C(individual)",
                data=df_model,
            ).fit()

            results.append(
                {
                    "gene": gene,
                    "variant": var_id,
                    "beta_genotype": model.params.get("genotype", np.nan),
                    "pval_genotype": model.pvalues.get("genotype", np.nan),
                    "beta_genox_med": model.params.get("genotype:C(activity)[T.medium]", np.nan),
                    "pval_genox_med": model.pvalues.get("genotype:C(activity)[T.medium]", np.nan),
                    "beta_genox_high": model.params.get("genotype:C(activity)[T.high]", np.nan),
                    "pval_genox_high": model.pvalues.get("genotype:C(activity)[T.high]", np.nan),
                    "beta_med": model.params.get("C(activity)[T.medium]", np.nan),
                    "pval_med": model.pvalues.get("C(activity)[T.medium]", np.nan),
                    "beta_high": model.params.get("C(activity)[T.high]", np.nan),
                    "pval_high": model.pvalues.get("C(activity)[T.high]", np.nan),
                }
            )

        except Exception as e:
            print(f"Model failed for {gene} - {var_id}: {e}")
            continue

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
@click.option('--job-cpu', help='Number of CPUs of Hail batch job', default=1)
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
            continue
        j = b.new_python_job(
            name=f'{cell_type}: {chromosome} {gene}',
        )
        j.cpu(job_cpu)
        j.memory(job_memory)
        j.storage(job_storage)
        j.call(process_gene, pheno_cov_dir, gene, chromosome, cell_type, pathway)

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter