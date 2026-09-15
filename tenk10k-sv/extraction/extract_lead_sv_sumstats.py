#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
Export the (lead SV x eGene x cell type) effect-size matrix (deliverable 3).

One row per lead SV, eGene and cell type -- INCLUDING the cell types where the association is not
significant. That inclusion is the whole point of the file: cell-type specificity, the cross-cell-
type correlation heatmap and any mash run are all statements about where an effect is ABSENT, and
absence cannot be read off a table that has been filtered to FDR < 5%.

Input is the significant-hit table (the FDR-filtered gene-level results, e.g. esvs.csv). That gives
the (federate_id, gene, chr) triples to look up; the betas then come from the full meta-analysis
output, unfiltered, across every cell type.

WHY THIS FILE IS NEEDED
-----------------------
run_gene_level_pvals_sv.py and run_gene_level_pvals_meta_sv.py both carry coeff/se for the single
lowest-raw-p locus per gene. raw_pval_extractor_sv.py keeps p-values only. So nothing currently
written out holds an effect size for a given locus across cell types, and there is no way to say
"this SV moves FCGR3A in NK but not in B cells" from the existing outputs -- the non-significant
cells are simply not represented.

Restricting to lead SVs keeps this small. ~5k lead SVs x their eGenes x 28 cell types is a few
hundred thousand rows, versus the hundreds of millions an unrestricted export would produce.

PRIVACY
-------
Every column is a regression summary statistic or a sample COUNT. No per-sample value is read or
written. `n_samples_tested_{1,2}` are per-cohort counts, already present in the meta output.

analysis-runner --dataset "tenk10k" --description "export lead SV sumstats" --access-level "test" \
    --output-dir "tenk10k-sv/sv/export/lead_sv_sumstats/v1" \
    python3 extract_lead_sv_sumstats.py \
    --lead-svs=gs://cpg-tenk10k-test/sv/export/esvs.csv \
    --meta-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/10pc/v1/meta_results \
    --cell-types=CD4_TCM,CD4_Naive,NK,CD8_TEM \
    --cis-window-size=1000000
"""

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import output_path

# Columns copied straight from the meta-analysis output. `gene` is renamed to `gene_name` on the way
# out so the export joins directly against the gene-level tables, which use `gene_name`.
PASSTHROUGH_COLUMNS = [
    'federate_id',
    'pair_class',
    'chr',
    'bioheart_vid',
    'bioheart_pos',
    'bioheart_svtype',
    'tob_vid',
    'tob_pos',
    'tob_svtype',
    'n_samples_tested_1',
    'n_samples_tested_2',
    'coeff_1',
    'se_1',
    'pval_1',
    'coeff_2',
    'se_2',
    'pval_2',
    'coeff_meta_fixed',
    'se_meta_fixed',
    'pval_meta_fixed',
    'coeff_meta_random',
    'se_meta_random',
    'pval_meta_random',
    'pval_q_meta',
    'i2_meta',
    'direction_concordant',
    'abs_alt_af_diff',
]

OUTPUT_COLUMNS = ['cell_type', 'gene_name', *PASSTHROUGH_COLUMNS]


def read_lead_set(lead_svs_path):
    """Significant-hit table -> {(chrom, gene): {federate_id, ...}}.

    Accepts either the raw esvs.csv (which uses `gene_name`) or any table carrying
    federate_id / gene / chr.
    """
    with to_path(lead_svs_path).open() as handle:
        lead = pd.read_csv(handle)

    gene_col = 'gene_name' if 'gene_name' in lead.columns else 'gene'
    missing = [c for c in ('federate_id', 'chr', gene_col) if c not in lead.columns]
    if missing:
        raise ValueError(
            f'{lead_svs_path} is missing column(s): {missing}. '
            'Expected federate_id, chr, and gene_name (or gene).',
        )

    wanted = defaultdict(set)
    for chrom, gene, fed in zip(lead['chr'], lead[gene_col], lead['federate_id']):
        chrom = chrom if str(chrom).startswith('chr') else f'chr{chrom}'
        wanted[(chrom, gene)].add(fed)

    print(
        f'{len(lead):,} significant rows -> {len(wanted):,} (chrom, gene) lookups, '
        f'{lead.federate_id.nunique():,} distinct lead SVs',
    )
    return dict(wanted)


def read_one(args):
    """Read one gene's meta results and keep the rows for its lead SVs.

    Returns (rows, status) where status is 'ok', 'missing' or 'no_match'. A gene that is significant
    in one cell type will legitimately have no meta result file in a cell type where it was never
    tested (too few expressing cells), so 'missing' is expected and is counted, not raised.
    """
    meta_dir, cell_type, chrom, gene, fed_ids, cis_window = args
    path = f'{meta_dir}/{cell_type}/{chrom}/{gene}_{cis_window}bp_meta_results.tsv'
    try:
        with to_path(path).open() as handle:
            res = pd.read_csv(handle, sep='\t')
    except (FileNotFoundError, OSError):
        return [], 'missing'
    if res.empty:
        return [], 'missing'

    res = res[res['federate_id'].isin(fed_ids)]
    if res.empty:
        return [], 'no_match'

    # A column absent from an older meta run is emitted as NA rather than dropped, so the schema is
    # stable across runs and the output can be concatenated.
    for col in PASSTHROUGH_COLUMNS:
        if col not in res.columns:
            res[col] = pd.NA

    res = res.assign(cell_type=cell_type, gene_name=gene)
    return res[OUTPUT_COLUMNS].values.tolist(), 'ok'


@click.command()
@click.option('--lead-svs', required=True, help='GCS path to the FDR-filtered hit table (e.g. esvs.csv)')
@click.option('--meta-dir', required=True, help='GCS path to meta_results from meta_runner_sv.py')
@click.option('--cell-types', required=True, help='Comma-separated cell types (give ALL of them)')
@click.option('--cis-window-size', default=1000000, type=int)
@click.option('--threads', default=32, type=int, help='concurrent GCS reads')
def main(lead_svs, meta_dir, cell_types, cis_window_size, threads):
    """Export lead-SV effect sizes across all cell types, significant or not.

    One output file per cell type, so a run can be sharded across several invocations -- the full
    28-cell-type sweep is ~7.7k GCS reads per cell type, which is more than one driver job wants to
    do serially. Concatenate the per-cell-type files afterwards.
    """
    wanted = read_lead_set(lead_svs)
    cell_type_list = [c.strip() for c in cell_types.split(',')]

    totals = {'ok': 0, 'missing': 0, 'no_match': 0}
    grand_rows = 0

    for cell_type in cell_type_list:
        jobs = [
            (meta_dir, cell_type, chrom, gene, fed_ids, cis_window_size)
            for (chrom, gene), fed_ids in wanted.items()
        ]
        rows_all = []
        with ThreadPoolExecutor(max_workers=threads) as executor:
            for rows, status in executor.map(read_one, jobs):
                totals[status] += 1
                rows_all.extend(rows)

        if not rows_all:
            print(f'  {cell_type}: no rows -- gene set never tested here? skipping file')
            continue

        out = pd.DataFrame(rows_all, columns=OUTPUT_COLUMNS)
        local = f'{cell_type}_lead_sv_sumstats.tsv.gz'
        # written locally then uploaded, matching sv_vcf_for_associatr.py -- pandas' own compression
        # handling against a remote file object is version-dependent, uploading a finished file is not
        out.to_csv(local, sep='\t', index=False, compression='gzip')
        gcs_output = output_path(f'lead_sv_sumstats/{local}', 'analysis')
        to_path(gcs_output).upload_from(local)
        grand_rows += len(out)
        print(f'  {cell_type}: {len(out):,} rows '
              f'({out.federate_id.nunique():,} SVs, {out.gene_name.nunique():,} genes) -> {gcs_output}')

    print(f'\ntotal {grand_rows:,} rows across {len(cell_type_list)} cell type(s)')
    print(f'gene lookups: {totals["ok"]:,} ok, {totals["missing"]:,} no result file '
          '(gene not tested in that cell type -- expected), '
          f'{totals["no_match"]:,} file present but lead SV absent')
    print('All columns are regression summary statistics or sample counts.')


if __name__ == '__main__':
    main()
