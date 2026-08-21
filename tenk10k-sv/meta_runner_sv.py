#!/usr/bin/env python3

"""
Meta-analyse associaTR SV eQTL summary statistics across the BioHEART and TOB cohorts.

Pooled effect sizes come from R's `meta` package (metagen) running in the `r-meta` image via rpy2 --
the same call, the same arguments and the same object fields as
str/associatr/meta_analysis/meta_runner.py -- so the meta columns are directly comparable with the
STR meta-analysis outputs.

The SV-specific part is the join. Unlike STRs, SV loci cannot be merged on (chrom, pos, motif):
the two cohorts were called independently, so the same underlying SV carries a different VARID and
often a different POS in each cohort. The join therefore goes through the FederateSV pairing table
produced by federate_overlap.ipynb:

    gs://cpg-tenk10k-test/sv/input_files/meta_analysis_set.csv

One row per federated variant, with `bioheart_vids` / `tob_vids` giving the per-cohort VARIDs and
`pair_class` recording how the two cohorts' SVTYPEs agreed. Two details make it work:

1. associaTR stores the VARID inside the `motif` column, as `<REF>-<ALT>-END:<end>-<VARID>`, e.g.
       N-<INS>-END:168956831-ALL_BATCHES.CHR1.FINAL_CLEANUP_INS_CHR1_2102
   so the VARID is whatever follows the `-END:<digits>-` segment.

2. That embedded VARID is UPPERCASED, while the VCF (and hence meta_analysis_set.csv) keeps the
   original mixed case (`all_batches.chr1.final_cleanup_INS_chr1_2102`). Both sides are upper-cased
   before matching.

Which pair classes are pooled
-----------------------------
`concordant_*` (both cohorts called the same SVTYPE) and `dup_cnv` (biallelic DUP in one cohort,
multiallelic CNV in the other) are pooled. For dup_cnv both codings measure effect per additional
copy -- CN = 2 + (dup allele count), and the +2 offset is absorbed by each cohort's own intercept --
so the betas share a scale and a direction. The biallelic coding saturates at CN 4, which attenuates
that cohort's beta toward zero; this costs power and understates effect magnitude but cannot create
false positives, so it is reported (via `pair_class`) rather than excluded. Use --exclude-dup-cnv to
run without them as a sensitivity check.

`del_cnv` pairs are ALWAYS dropped: dosage_DEL = 2 - CN runs opposite to dosage_CNV, so pooling
would cancel real signal. They should already be absent from meta_analysis_set.csv; if any are
present the script warns loudly and drops them.

Output: one TSV per gene, mirroring the STR meta-analysis layout.
    {output-dir}/meta_results/{cell_type}/{chrom}/{gene}_{cis_window}bp_meta_results.tsv

analysis-runner --dataset "tenk10k" --description "SV eQTL meta-analysis" --access-level "test" \
    --output-dir "tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/10pc/v1" \
    --image 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:189f33cce2381aa8623acd2f39006455deae2a6f-hail-0928b141d1dd3ca4de530d50943d2673305ea7bb' \
    python3 meta_runner_sv.py \
    --results-dir-1=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/sensitivity_analysis/bioheart_n968/9pc/results/v1 \
    --results-dir-2=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/sensitivity_analysis/tob_n935/10pc/results/v1 \
    --meta-set=gs://cpg-tenk10k-test/sv/input_files/meta_analysis_set.csv \
    --cell-types=CD8_TEM \
    --chromosomes=chr1

"""

import json
import re

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path

# associaTR RU/motif field is '<REF>-<ALT>-END:<end>-<VARID>'; the VARID is the tail.
_VARID_RE = re.compile(r'-END:\d+-(.+)$')

_POOLABLE_PREFIX = 'concordant_'
_POOLABLE_EXTRA = 'dup_cnv'

# metagen object field -> output column name. Field names match meta_runner.py's m.gen$... accessors
# so the pooled columns line up with the STR meta-analysis outputs. `meta` >= 5 renamed the
# fixed-effect fields from *.fixed to *.common and kept the old names as aliases; both are tried.
_METAGEN_FIELDS = (
    (('TE.random',), 'coeff_meta_random'),
    (('seTE.random',), 'se_meta_random'),
    (('pval.Q',), 'pval_q_meta'),
    (('Q',), 'q_meta'),
    (('tau',), 'tau_meta'),
    (('I2',), 'i2_meta'),
    (('H',), 'h_meta'),
    (('pval.random',), 'pval_meta_random'),
    (('pval.fixed', 'pval.common'), 'pval_meta_fixed'),
    (('TE.fixed', 'TE.common'), 'coeff_meta_fixed'),
    (('seTE.fixed', 'seTE.common'), 'se_meta_fixed'),
    (('lower.random',), 'lowerCI_meta_random'),
    (('upper.random',), 'upperCI_meta_random'),
)


def _varid_from_motif(motif):
    """'N-<DUP>-END:168957979-ALL_BATCHES...DUP_CHR1_912' -> 'ALL_BATCHES...DUP_CHR1_912' (upper).

    Returns None when the motif does not carry a recognisable VARID.
    """
    match = _VARID_RE.search(str(motif))
    return match.group(1).strip().upper() if match else None


def _parse_af(af_json):
    """allele_frequency JSON -> (alt_af, mean_dosage).

    `alt_af` is only defined for a biallelic 0/1 dosage coding; multiallelic CNV records get NaN.
    `mean_dosage` = sum(dosage * freq) is defined for any coding, but it differs by ~2 between a
    biallelic DUP coding and a CN coding at the same locus, so only compare within a pair class.
    """
    nan = float('nan')
    try:
        freqs = {float(k): float(v) for k, v in json.loads(af_json).items()}
    except (TypeError, ValueError, json.JSONDecodeError):
        return nan, nan
    mean_dosage = sum(dosage * freq for dosage, freq in freqs.items())
    alt_af = freqs.get(1.0, nan) if set(freqs) <= {0.0, 1.0} else nan
    return alt_af, mean_dosage


def _read_summary_stats(path):
    """Read one associaTR gene TSV, keeping tested loci keyed on the upper-cased VARID.

    Returns None if the file is missing/unreadable, so a gene present in only one cohort is skipped
    rather than failing the whole job.
    """
    try:
        with to_path(path).open() as handle:
            df = pd.read_csv(handle, sep='\t')
    except (FileNotFoundError, OSError):
        return None
    if df.empty:
        return None

    # column names embed cell type and gene (e.g. coeff_CD8_TEM_chr1_ENSG00000000457), so locate
    # them by prefix rather than reconstructing the name
    renames = {}
    for prefix, canonical in (('coeff_', 'coeff'), ('se_', 'se'), ('p_', 'pval')):
        hits = [c for c in df.columns if c.startswith(prefix)]
        if len(hits) != 1:
            raise ValueError(f'expected exactly one {prefix}* column in {path}, found {hits}')
        renames[hits[0]] = canonical
    df = df.rename(columns=renames)

    # associaTR flags loci it could not test
    df = df[df['locus_filtered'].astype(str) == 'False']
    df['varid'] = df['motif'].map(_varid_from_motif)
    df = df[df['varid'].notna()]

    # zero/missing SEs would blow up metagen's inverse-variance weights
    df = df[df['se'].notna() & (df['se'] > 0) & df['coeff'].notna()]

    af = df['allele_frequency'].map(_parse_af)
    df['alt_af'] = [a for a, _ in af]
    df['mean_dosage'] = [m for _, m in af]

    return df.drop_duplicates(subset='varid', keep=False).set_index('varid')


def _load_pairs(meta_set_path, include_dup_cnv=True):
    """meta_analysis_set.csv -> DataFrame of poolable pairs, keyed by upper-cased cohort VARIDs."""
    with to_path(meta_set_path).open() as handle:
        pairs = pd.read_csv(handle)

    n_del_cnv = int((pairs['pair_class'] == 'del_cnv').sum())
    if n_del_cnv:
        print(
            f'!! WARNING: dropping {n_del_cnv} del_cnv pair(s) -- dosage_DEL = 2 - CN runs opposite '
            f'to dosage_CNV, so pooling them would cancel real signal. Was {meta_set_path} the '
            f'filtered meta-analysis set?',
        )

    keep = pairs['pair_class'].str.startswith(_POOLABLE_PREFIX)
    if include_dup_cnv:
        keep = keep | (pairs['pair_class'] == _POOLABLE_EXTRA)
    pairs = pairs[keep].copy()

    # every row of the current table is 1:1, but tolerate comma-joined multi-call rows by taking
    # the first VARID and reporting how many rows that simplification touched
    for col in ('bioheart_vids', 'tob_vids'):
        multi = pairs[col].astype(str).str.contains(',')
        if multi.any():
            print(f'  note: {int(multi.sum())} pair(s) have multiple {col}; using the first of each')
        pairs[col] = pairs[col].astype(str).str.split(',').str[0].str.strip()

    pairs['bh_key'] = pairs['bioheart_vids'].str.upper()
    pairs['tob_key'] = pairs['tob_vids'].str.upper()
    return pairs


def run_meta_for_chrom(
    results_dir_1,
    results_dir_2,
    meta_set_path,
    cell_type,
    chromosome,
    cis_window,
    include_dup_cnv,
):
    """Meta-analyse every gene tested in both cohorts for one cell type / chromosome."""
    import rpy2.robjects as ro

    from cpg_utils.hail_batch import output_path

    # load once per job rather than once per gene
    ro.r('library(meta)')

    def metagen_pair(coeff_1, se_1, coeff_2, se_2):
        """One metagen() call for a two-study meta-analysis; returns the STR-named fields.

        Identical call and arguments to meta_runner.py:
            metagen(result_df$coeff, result_df$se, random = TRUE, method.tau = "DL")
        """
        ro.globalenv['te'] = ro.FloatVector([coeff_1, coeff_2])
        ro.globalenv['sete'] = ro.FloatVector([se_1, se_2])
        m_gen = ro.r('metagen(te, sete, random = TRUE, method.tau = "DL")')
        available = set(m_gen.names)
        out = {}
        for candidates, column in _METAGEN_FIELDS:
            field = next((c for c in candidates if c in available), None)
            value = m_gen.rx2(field) if field else None
            out[column] = float(value[0]) if value is not None and len(value) else float('nan')
        return out

    pairs = _load_pairs(meta_set_path, include_dup_cnv)
    # restrict to variants on this chromosome to keep the per-gene joins small
    pairs = pairs[pairs['federate_chrom'] == chromosome]
    print(f'{cell_type} {chromosome}: {len(pairs):,} poolable federated variants')
    if pairs.empty:
        return

    suffix = f'_{cis_window}bp.tsv'
    genes_1 = {p.name[: -len(suffix)] for p in to_path(f'{results_dir_1}/{cell_type}/{chromosome}').glob(f'*{suffix}')}
    genes_2 = {p.name[: -len(suffix)] for p in to_path(f'{results_dir_2}/{cell_type}/{chromosome}').glob(f'*{suffix}')}
    genes = sorted(genes_1 & genes_2)
    print(f'  {len(genes):,} genes tested in both cohorts ({len(genes_1):,} / {len(genes_2):,} per cohort)')

    n_written = n_skipped = 0
    # genes that yield no output file, tracked so written + skipped + these == len(genes)
    unreadable = []  # a cohort TSV was missing, empty, or had no testable loci left
    no_pairs = []  # both cohorts had loci, but no federated pair had both sides poolable
    for gene in genes:
        out_path = output_path(
            f'meta_results/{cell_type}/{chromosome}/{gene}_{cis_window}bp_meta_results.tsv',
            'analysis',
        )
        if to_path(out_path).exists():
            n_skipped += 1
            continue

        d1 = _read_summary_stats(f'{results_dir_1}/{cell_type}/{chromosome}/{gene}{suffix}')
        d2 = _read_summary_stats(f'{results_dir_2}/{cell_type}/{chromosome}/{gene}{suffix}')
        if d1 is None or d2 is None:
            unreadable.append(gene)
            continue

        # join each cohort's summary stats onto the federated pairing table by upper-cased VARID
        merged = pairs[pairs['bh_key'].isin(d1.index) & pairs['tob_key'].isin(d2.index)]
        if merged.empty:
            no_pairs.append(gene)
            continue

        rows = []
        for pair in merged.itertuples(index=False):
            r1, r2 = d1.loc[pair.bh_key], d2.loc[pair.tob_key]
            row = {
                'federate_id': pair.federate_id,
                'pair_class': pair.pair_class,
                'cell_type': cell_type,
                'chr': pair.federate_chrom,
                'gene': gene,
                'bioheart_vid': pair.bioheart_vids,
                'bioheart_pos': pair.bioheart_pos,
                'bioheart_svtype': pair.bioheart_svtype,
                'tob_vid': pair.tob_vids,
                'tob_pos': pair.tob_pos,
                'tob_svtype': pair.tob_svtype,
                'n_samples_tested_1': r1['n_samples_tested'],
                'n_samples_tested_2': r2['n_samples_tested'],
                'coeff_1': r1['coeff'],
                'se_1': r1['se'],
                'pval_1': r1['pval'],
                'r2_1': r1['regression_R^2'],
                'motif_1': r1['motif'],
                'alleles_1': r1['alleles'],
                'allele_frequency_1': r1['allele_frequency'],
                'alt_af_1': r1['alt_af'],
                'coeff_2': r2['coeff'],
                'se_2': r2['se'],
                'pval_2': r2['pval'],
                'r2_2': r2['regression_R^2'],
                'motif_2': r2['motif'],
                'alleles_2': r2['alleles'],
                'allele_frequency_2': r2['allele_frequency'],
                'alt_af_2': r2['alt_af'],
            }
            row.update(metagen_pair(r1['coeff'], r1['se'], r2['coeff'], r2['se']))
            row['direction_concordant'] = (r1['coeff'] > 0) == (r2['coeff'] > 0)
            # a large AF gap between two variants FederateSV called the same is a merge-quality flag,
            # not something to correct for; only defined when both cohorts are biallelic
            row['abs_alt_af_diff'] = abs(r1['alt_af'] - r2['alt_af'])
            rows.append(row)

        out_df = pd.DataFrame(rows).sort_values('pval_meta_fixed')
        with to_path(out_path).open('w') as handle:
            out_df.to_csv(handle, sep='\t', index=False)
        n_written += 1

    print(f'  wrote {n_written:,} gene result files ({n_skipped:,} already present, skipped)')
    if no_pairs:
        shown = ', '.join(no_pairs[:10]) + (' ...' if len(no_pairs) > 10 else '')
        print(f'  {len(no_pairs):,} gene(s) had no federated pair in the cis window: {shown}')
    if unreadable:
        shown = ', '.join(unreadable[:10]) + (' ...' if len(unreadable) > 10 else '')
        print(f'  {len(unreadable):,} gene(s) had an empty/unreadable cohort TSV: {shown}')
    accounted = n_written + n_skipped + len(no_pairs) + len(unreadable)
    if accounted != len(genes):
        print(f'  !! accounting gap: {len(genes):,} genes seen but {accounted:,} accounted for')


@click.command()
@click.option('--results-dir-1', required=True, help='associaTR results dir for cohort 1 (BioHEART)')
@click.option('--results-dir-2', required=True, help='associaTR results dir for cohort 2 (TOB)')
@click.option(
    '--meta-set',
    default='gs://cpg-tenk10k-test/sv/input_files/meta_analysis_set.csv',
    help='FederateSV pairing table (output of federate_overlap.ipynb)',
)
@click.option('--cell-types', required=True, help='comma-separated cell types')
@click.option('--chromosomes', required=True, help='comma-separated chromosomes')
@click.option('--cis-window', default=1000000, type=int, help='cis window in the result filenames')
@click.option(
    '--exclude-dup-cnv',
    is_flag=True,
    help='drop the DUP-vs-CNV pairs as well (sensitivity check on the biallelic/multiallelic coding)',
)
@click.option('--job-cpu', default=1, type=float, help='CPUs per job')
@click.option('--max-parallel-jobs', default=500, type=int, help='concurrency cap')
@click.option('--always-run', is_flag=True, help='set jobs to always run')
def main(
    results_dir_1,
    results_dir_2,
    meta_set,
    cell_types,
    chromosomes,
    cis_window,
    exclude_dup_cnv,
    job_cpu,
    max_parallel_jobs,
    always_run,
):
    """Fan out one job per cell type / chromosome; each loops over the genes tested in both cohorts."""
    from cpg_utils.hail_batch import get_batch, image_path

    include_dup_cnv = not exclude_dup_cnv

    # report what is being pooled before spending any compute
    pairs = _load_pairs(meta_set, include_dup_cnv)
    print(f'Pooling {len(pairs):,} federated variants by pair_class:')
    for cls, n in pairs['pair_class'].value_counts().items():
        print(f'  {cls}: {n:,}')
    if exclude_dup_cnv:
        print('  (dup_cnv excluded by --exclude-dup-cnv)')

    b = get_batch(name='SV eQTL meta-analysis')

    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """Limit how many jobs run at once."""
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    for cell_type in cell_types.split(','):
        for chromosome in chromosomes.split(','):
            j = b.new_python_job(name=f'meta_{cell_type}_{chromosome}')
            j.image(image_path('r-meta'))
            j.cpu(job_cpu)
            if always_run:
                j.always_run()
            j.call(
                run_meta_for_chrom,
                results_dir_1,
                results_dir_2,
                meta_set,
                cell_type,
                chromosome,
                cis_window,
                include_dup_cnv,
            )
            manage_concurrency_for_job(j)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
