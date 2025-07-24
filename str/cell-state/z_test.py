#!/usr/bin/env python3

"""
This script runs z-test to compare TR-eQTL effect sizes between activity bins.

analysis-runner --dataset "tenk10k" --description "z-test for context-spec eQTLs" --access-level "test" \
    --output-dir "str/cellstate/stratified/meta_results/GOBP_MULTI" z_test.py

"""

import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, init_batch, output_path

import click
import pandas as pd


def z_test_runner(meta_dir, pathway, cell_type, chromosome, pos, end,motif, gene, ref_activity_level):
    import pandas as pd
    import numpy as np
    from scipy.stats import norm
    """
    Run Z-tests to compare TR-eQTL effect sizes between activity bins.

    Parameters
    ----------
    meta_dir : str
        Root directory containing meta-analysis result files.
    pathway : str
        Pathway or activity score name.
    cell_type : str
        Major or fine cell type.
    chromosome : str
        Chromosome name (e.g. 'chr1').
    pos : int
        TR genomic position.
    motif : str
        TR motif string.
    gene : str
        Gene symbol.
    ref_activity_level : str
        Activity level to treat as the reference bin ('low', 'medium', or 'high').

    Returns
    -------
    results : list of dict
        Each dict contains the comparison name, delta_beta, Z, P, and metadata.
    """
    activity_levels = ['low', 'medium', 'high']
    assert ref_activity_level in activity_levels, "Invalid reference activity level"

    other_levels = [lvl for lvl in activity_levels if lvl != ref_activity_level]

    def load_and_filter(activity_level):
        filepath = f"{meta_dir}/{pathway}/{activity_level}/{cell_type}/{chromosome}/{gene}_100000bp_meta_results.tsv"
        df = pd.read_csv(filepath, sep='\t')
        df['motif_len'] = df['motif'].str.len()
        df['end'] = (df['pos'].astype(float) + df['ref_len'].astype(float) * df['motif_len'].astype(float)).round().astype(int)
        row = df[(df['pos'] == pos) & (df['motif'] == motif) &(df['end'] == end)]
        return row

    ref_row = load_and_filter(ref_activity_level)
    row_1 = load_and_filter(other_levels[0])
    row_2 = load_and_filter(other_levels[1])

    # Check for missing entries
    if ref_row.empty or row_1.empty or row_2.empty:
        return []  # No result if any activity level is missing for this variant

    ref_beta = ref_row['coeff_meta_fixed'].values[0]
    ref_se = ref_row['se_meta_fixed'].values[0]

    results = []

    for other_level, other_row in zip(other_levels, [row_1, row_2]):
        other_beta = other_row['coeff_meta_fixed'].values[0]
        other_se = other_row['se_meta_fixed'].values[0]

        delta_beta = other_beta - ref_beta
        denom = np.sqrt(ref_se**2 + other_se**2)

        if denom == 0:
            continue  # Avoid division by zero

        z = delta_beta / denom
        pval = 2 * norm.sf(abs(z))

        results.append({
            "chrom": chromosome,
            "pos": pos,
            "end": end,
            "motif": motif,
            "gene": gene,
            "cell_type": cell_type,
            "pathway": pathway,
            "comparison": f"{other_level}_vs_{ref_activity_level}",
            "beta_diff": delta_beta,
            "z_score": z,
            "p_value": pval,
            "ref_beta": ref_beta,
            "ref_se": ref_se,
            "other_beta": other_beta,
            "other_se": other_se
        })
    # Comparison between the two other  non-ref bins
    level_a, level_b = other_levels
    row_a, row_b = row_1, row_2

    beta_a = row_a['coeff_meta_fixed'].values[0]
    se_a = row_a['se_meta_fixed'].values[0]
    beta_b = row_b['coeff_meta_fixed'].values[0]
    se_b = row_b['se_meta_fixed'].values[0]

    delta_beta_ab = beta_a - beta_b
    denom_ab = np.sqrt(se_a**2 + se_b**2)

    if denom_ab != 0:
        z_ab = delta_beta_ab / denom_ab
        pval_ab = 2 * norm.sf(abs(z_ab))

        results.append({
            "chrom": chromosome,
            "pos": pos,
            "end": end,
            "motif": motif,
            "gene": gene,
            "cell_type": cell_type,
            "pathway": pathway,
            "comparison": f"{level_a}_vs_{level_b}",
            "beta_diff": delta_beta_ab,
            "z_score": z_ab,
            "p_value": pval_ab,
            "ref_beta": beta_b,  # treat 'b' as reference here
            "ref_se": se_b,
            "other_beta": beta_a,
            "other_se": se_a
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv(
        f"{meta_dir}/{pathway}/z_test/{cell_type}/{chromosome}/{gene}_z_test.tsv",
        index=False,
        sep='\t')

@click.option('--meta-dir', help='Directory containing meta-analysis results', default= 'gs://cpg-tenk10k-test-analysis/str/cellstate/stratified/meta_results')
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

        if to_path(f'{meta_dir}/{pathway}/z_test/{cell_type}/{chromosome}/{gene}_z_test.tsv').exists():
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
