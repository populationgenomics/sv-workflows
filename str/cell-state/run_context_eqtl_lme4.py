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
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from cpg_utils.hail_batch import output_path

    # Load phenotype + covariate data
    df = pd.read_csv(f'{pheno_cov_dir}/{pathway}/{cell_type}/{chromosome}/{gene}_pheno_cov.csv')

    # Load TR dosage matrix
    dosage = pd.read_csv(f'gs://cpg-tenk10k-test/str/cellstate/input_files/dosages/{chromosome}/{gene}_dosages.csv')
    dosage['sample'] = 'CPG' + dosage['sample'].astype(str)

    # Merge by sample ID
    df = df.merge(dosage, left_on='sample_id', right_on='sample', how='inner')
    df['activity'] = pd.Categorical(df['activity'], categories=['low', 'medium', 'high'], ordered=False)
    df['individual'] = df['sample_id'].astype(str)  # for grouping in lme4

    variant_cols = [col for col in dosage.columns if col != 'sample']
    results = []


    # Setup R libraries
    ro.r('library(lme4)')
    ro.r('library(broom.mixed)')

    for var_id in variant_cols:
        df['genotype'] = df[var_id]
        df_model = df.dropna(subset=["gene_inverse_normal", "genotype"])

        with (ro.default_converter + pandas2ri.converter).context():
            df_r = ro.conversion.py2rpy(df_model)
            ro.globalenv['df'] = df_r
            ro.globalenv['var'] = var_id
            ro.globalenv['gene'] = gene

            try:
                ro.r(
                    """
                    df$activity <- factor(df$activity, levels=c('low','medium','high'))
                    df$individual <- factor(df$individual)

                    model <- lmer(gene_inverse_normal ~ genotype * activity +
                                  sex + age +
                                  score_1 + score_2 + score_3 + score_4 + score_5 + score_6 +
                                  score_7 + score_8 + score_9 + score_10 + score_11 + score_12 +
                                  rna_PC1 + rna_PC2 + rna_PC3 + rna_PC4 + rna_PC5 + rna_PC6 +
                                  (1 | individual),
                                  data = df, REML = FALSE)

                    out <- broom.mixed::tidy(model, effects = "fixed")
                    out$variant <- var
                    out$gene <- gene
                    """,
                )
                with (ro.default_converter + pandas2ri.converter).context():
                    out_df = ro.conversion.rpy2py(ro.r('out'))
                    results.append(out_df)
            except Exception as e:
                print(f"Model failed for {gene} - {var_id}: {e}")
                continue

    if results:
        results_df = pd.concat(results)
        results_df['cell_type'] = cell_type
        results_df['chromosome'] = chromosome

        # Write to GCS
        results_df.to_csv(
            output_path(f'context_eqtl_results/{pathway}/{cell_type}/{chromosome}/{gene}_context_eqtl_results.csv'),
            index=False,
        )

@click.option('--input-gene-list-dir', default='gs://cpg-tenk10k-test/str/cellstate/input_files/tob/scRNA_gene_lists')
@click.option('--pheno-cov-dir', default='gs://cpg-tenk10k-test/str/cellstate/input_files/tob/pheno_cov_csv')
@click.option('--cell-type', default='B_naive', help='Cell type to process')
@click.option('--pathway', default='GOBP_MULTI_MULTICELLULAR_ORGANISM_PROCESS_subtype', help='Pathway to process')
@click.option('--chromosome', default='chr22', help='Chromosome to process')
@click.option('--job-storage', help='Storage of the batch job eg 30G', default='4G')
@click.option('--job-memory', help='Memory of the batch job', default='standard')
@click.option('--job-cpu', help='Number of CPUs of Hail batch job', default=0.25)
@click.command()
def main(
    input_gene_list_dir, pheno_cov_dir, cell_type, pathway, chromosome, job_storage, job_memory, job_cpu
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