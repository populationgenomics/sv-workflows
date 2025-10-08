#!/usr/bin/env python3

"""
This script runs z-test to compare TR-eQTL effect sizes between activity bins.

analysis-runner --dataset "tenk10k" --description "z-test for context-spec eQTLs" --access-level "test" \
    --output-dir "str/cellstate/stratified/meta_results/GOBP_MULTI" z_test.py \
    --eqtls-to-test=gs://cpg-tenk10k-test-analysis/str/cellstate/meanpool/stratified/meta_results/GOBP_MULTI_MULTICELLULAR_ORGANISM_PROCESS_subtype/eqtls_to_test.csv

"""

import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path

import click
import pandas as pd


def z_test_runner(meta_dir, pathway, cell_type, chromosome, pos, end, motif, gene, ref_activity_level):
    import pandas as pd
    import numpy as np
    from scipy.stats import norm
    import os

    activity_levels = ['low', 'medium', 'high']
    assert ref_activity_level in activity_levels, "Invalid reference activity level"
    other_levels = [lvl for lvl in activity_levels if lvl != ref_activity_level]

    def load_and_filter(activity_level):
        filepath = f"{meta_dir}/{pathway}/{activity_level}/{cell_type}/{chromosome}/{gene}_100000bp_meta_results.tsv"
        df = pd.read_csv(filepath, sep='\t')
        df['motif_len'] = df['motif'].str.len()
        df['end'] = (df['pos'].astype(float) + df['ref_len'].astype(float) * df['motif_len'].astype(float)).round().astype(int)
        row = df[(df['pos'] == pos) & (df['motif'] == motif) & (df['end'] == end)]
        return row

    rows = {}
    for level in activity_levels:
        row = load_and_filter(level)
        if row.empty:
            return []  # Skip this variant if missing in any bin
        rows[level] = row

    # Extract values
    beta_low = rows['low']['coeff_meta_fixed'].values[0]
    beta_medium = rows['medium']['coeff_meta_fixed'].values[0]
    beta_high = rows['high']['coeff_meta_fixed'].values[0]
    se_low = rows['low']['se_meta_fixed'].values[0]
    se_medium = rows['medium']['se_meta_fixed'].values[0]
    se_high = rows['high']['se_meta_fixed'].values[0]

    # Prepare results for z-tests
    results = []

    def compare(beta1, se1, beta2, se2, label1, label2, ref_beta_for_row):
        delta = beta1 - beta2
        denom = np.sqrt(se1 ** 2 + se2 ** 2)
        if denom == 0:
            return None
        z = delta / denom
        pval = 2 * norm.sf(abs(z))
        return {
            "chrom": chromosome,
            "pos": pos,
            "end": end,
            "motif": motif,
            "gene": gene,
            "cell_type": cell_type,
            "pathway": pathway,
            "comparison": f"{label1}_vs_{label2}",
            "beta_diff": delta,
            "z_score": z,
            "p_value": pval,
            "ref_beta": ref_beta_for_row,
            "other_beta": beta1 if ref_beta_for_row == beta2 else beta2,
            "ref_se": se1 if ref_beta_for_row == beta2 else se2,
            "other_se": se1 if ref_beta_for_row == beta1 else se2
        }

    # Run pairwise comparisons
    comp1 = compare(rows[other_levels[0]]['coeff_meta_fixed'].values[0],
                    rows[other_levels[0]]['se_meta_fixed'].values[0],
                    rows[ref_activity_level]['coeff_meta_fixed'].values[0],
                    rows[ref_activity_level]['se_meta_fixed'].values[0],
                    other_levels[0], ref_activity_level,
                    rows[ref_activity_level]['coeff_meta_fixed'].values[0])
    if comp1:
        results.append(comp1)

    comp2 = compare(rows[other_levels[1]]['coeff_meta_fixed'].values[0],
                    rows[other_levels[1]]['se_meta_fixed'].values[0],
                    rows[ref_activity_level]['coeff_meta_fixed'].values[0],
                    rows[ref_activity_level]['se_meta_fixed'].values[0],
                    other_levels[1], ref_activity_level,
                    rows[ref_activity_level]['coeff_meta_fixed'].values[0])
    if comp2:
        results.append(comp2)

    # Comparison between the two non-ref bins
    comp3 = compare(rows[other_levels[0]]['coeff_meta_fixed'].values[0],
                    rows[other_levels[0]]['se_meta_fixed'].values[0],
                    rows[other_levels[1]]['coeff_meta_fixed'].values[0],
                    rows[other_levels[1]]['se_meta_fixed'].values[0],
                    other_levels[0], other_levels[1],
                    rows[other_levels[1]]['coeff_meta_fixed'].values[0])
    if comp3:
        results.append(comp3)

    results_df = pd.DataFrame(results)

    out_dir = f"{meta_dir}/{pathway}/z_test/{cell_type}/{chromosome}"

    # Save z-test results
    results_df.to_csv(
        f"{out_dir}/{gene}_{pos}_{end}_{motif}_z_test.tsv",
        index=False,
        sep='\t'
    )

    # Save single-row summary of all effect sizes
    beta_se_df = pd.DataFrame([{
        "chrom": chromosome,
        "pos": pos,
        "end": end,
        "motif": motif,
        "gene": gene,
        "coeff_low": beta_low,
        "coeff_medium": beta_medium,
        "coeff_high": beta_high,
        "se_low": se_low,
        "se_medium": se_medium,
        "se_high": se_high
    }])
    beta_se_df.to_csv(
        f"{out_dir}/{gene}_{pos}_{end}_{motif}_effect_sizes.tsv",
        index=False,
        sep='\t'
    )

    return results


@click.option('--meta-dir', help='Directory containing meta-analysis results', default= 'gs://cpg-tenk10k-test-analysis/str/cellstate/meanpool/stratified/meta_results')
@click.option(
    '--pathway', help='GOBP Pathway name for the analysis', default='GOBP_MULTI_MULTICELLULAR_ORGANISM_PROCESS_subtype'
)
@click.option('--eqtls-to-test', help='File containing eQTLs to test', default='gs://cpg-tenk10k-test-analysis/str/cellstate/stratified/meta_results/GOBP_MULTI_MULTICELLULAR_ORGANISM_PROCESS_subtype/eqtls.csv')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs to run', default=1000)
@click.command()
def main(meta_dir,pathway,eqtls_to_test, max_parallel_jobs):
    """
    Run meta-analysis for eQTLs
    """
    b = get_batch(name='Z-test for context-spec eQTLs')
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    eqtls = pd.read_csv(eqtls_to_test)
    for row in eqtls.itertuples():
        cell_type = row.cell_type
        chromosome = row.chr
        pos = row.pos
        end = row.end
        motif = row.motif
        gene = row.gene_name
        ref_activity_level = row.activity_level

        if to_path(f'{meta_dir}/{pathway}/z_test/{cell_type}/{chromosome}/{gene}_{pos}_{end}_{motif}_z_test.tsv').exists():
            print(f"Z-test results for {gene} already exists, skipping.")
            #continue

        # Create a job for each eQTL
        j = b.new_python_job(
            name=f'Z-test for {cell_type} {chromosome} {pos} {motif} {gene}',
        )
        j.cpu(0.25)

        j.call(
            z_test_runner,
            meta_dir,
            pathway,
            cell_type,
            chromosome,
            pos,
            end,
            motif,
            gene,
            ref_activity_level,
        )

        manage_concurrency_for_job(j)
    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
