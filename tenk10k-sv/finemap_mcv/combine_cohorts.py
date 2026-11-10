#!/usr/bin/env python3
"""
Match SuSiE credible sets across the per-cohort (BioHEART, TOB) and pooled runs of susie_mcv_sv_runner.py.

SNPs and TRs carry the same IDs in every run (jointly called). Every SV is identified by its federated PAIR,
'FED:<BIOHEART_VARID>|<TOB_VARID>' (upper-cased, as meta_runner_sv.py matches them): per-cohort VARIDs are looked up
in the pairing table, pooled SVs carry bioheart_id / tob_id. federate_id is not used as the key because it is not
unique in the fixed pairing tables. An SV with no partner keeps a cohort-tagged ID and never matches across runs.

Two credible sets for the same gene x cell type match if they share at least one member after that mapping. Each
pair of runs -- bioheart/tob, pooled/bioheart, pooled/tob -- is matched one-to-one, largest overlap (Jaccard)
first.

Output, one row per credible set in any run:
    gene, cell_type, cohort, cs_id, cs_size, min_abs_corr, lead_variant, lead_type, alpha_SNP/TR/SV, best_sv,
    best_sv_pair (FED:<BH VARID>|<TOB VARID>), best_sv_alpha, best_sv_pip, best_sv_overlapping_trs,
    lead_esv_in_cs                       the gene x cell type's lead eSV (SV q-value file) is in this CS
    for each other run R in {bioheart, tob, pooled}:
        R_cs_id, R_jaccard               the matched CS in run R, if any
        R_sv_shared                      the same federated SV is a member of both credible sets
        R_sv_led                         ... and it leads both

Usage (reads gs:// paths directly):
    python3 combine_cohorts.py \
        --results-dir gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/finemap_mcv/bioheart_n968_and_tob_n935/7pc/v2 \
        --out finemap_cs_matched.tsv.gz

--meta-set and --qval-dir default to the v2 meta-analysis; they must match the run that was fine-mapped.
"""
import argparse
import io
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from cpg_utils import to_path

META_SET = 'gs://cpg-tenk10k-test/sv/input_files/meta_analysis_set_v2_fixed_F_filtered.csv'
QVAL_DIR = ('gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/'
            'v2_F_filtered/fdr_qvals/using_acat')
RUNS = ('bioheart', 'tob', 'pooled')


def read_tsv(p):
    with to_path(p).open('rb') as fh:
        raw = fh.read()
    if not raw:
        return pd.DataFrame()
    try:
        return pd.read_csv(io.BytesIO(raw), sep='\t', compression='gzip' if str(p).endswith('.gz') else None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def read_dir(d):
    files = list(to_path(d).glob('*.tsv.gz'))
    with ThreadPoolExecutor(16) as ex:
        dfs = [x for x in ex.map(read_tsv, files) if len(x)]
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def pair_key(bh, tob):
    return f'FED:{str(bh).split(",")[0].strip().upper()}|{str(tob).split(",")[0].strip().upper()}'


def federate_map(meta_set):
    """(cohort, upper-cased VARID) -> pair key, over the whole table (any pair class), first VARID per cell."""
    m = pd.read_csv(meta_set, usecols=['bioheart_vids', 'tob_vids']).dropna()
    out = {}
    for bh, tob in zip(m.bioheart_vids, m.tob_vids):
        k = pair_key(bh, tob)
        out[('bioheart', str(bh).split(',')[0].strip().upper())] = k
        out[('tob', str(tob).split(',')[0].strip().upper())] = k
    return out


def lead_esvs(qval_dir):
    """(gene, cell_type, run) -> set of lead eSV IDs: upper-cased per-cohort VARIDs, pair keys for pooled."""
    out = {}
    for f in to_path(qval_dir).glob('*_qval.tsv'):
        ct = f.name.removesuffix('_qval.tsv')
        d = pd.read_csv(f, sep='\t', usecols=['gene_name', 'qval', 'bioheart_vid', 'tob_vid'])
        for r in d[d.qval < 0.05].itertuples():
            bh, tob = str(r.bioheart_vid).split(','), str(r.tob_vid).split(',')
            out[(r.gene_name, ct, 'bioheart')] = {v.upper() for v in bh}
            out[(r.gene_name, ct, 'tob')] = {v.upper() for v in tob}
            out[(r.gene_name, ct, 'pooled')] = {pair_key(b, t) for b, t in zip(bh, tob)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results-dir', required=True)
    ap.add_argument('--meta-set', default=META_SET)
    ap.add_argument('--qval-dir', default=QVAL_DIR)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    cs = read_dir(f'{a.results_dir}/cs').reset_index(drop=True)
    var = read_dir(f'{a.results_dir}/variants')
    print(f'{len(cs):,} credible sets over {cs.groupby(["gene", "cell_type"]).ngroups:,} gene x cell type pairs; '
          f'runs: {cs.cohort.value_counts().to_dict()}')
    fmap = federate_map(a.meta_set)
    leads = lead_esvs(a.qval_dir)
    svs = var[var.variant_type == 'SV']
    sv_ids = set(zip(svs.cohort, svs.variant_id))
    pooled_key = {v: pair_key(b, t) for v, b, t in svs.loc[svs.cohort == 'pooled', ['variant_id', 'bioheart_id', 'tob_id']]
                  .drop_duplicates().itertuples(index=False)}

    def canon(v, run):
        if (run, v) not in sv_ids:
            return v
        if run == 'pooled':
            return pooled_key[v]
        return fmap.get((run, v.upper()), f'{run.upper()}:{v}')

    cs['members_c'] = [frozenset(canon(v, c) for v in m.split(',')) for m, c in zip(cs.members, cs.cohort)]
    cs['lead_c'] = [canon(v, c) for v, c in zip(cs.lead_variant, cs.cohort)]
    cs['best_sv_pair'] = [canon(v, c) if isinstance(v, str) else None for v, c in zip(cs.best_sv, cs.cohort)]
    cs['lead_esv_in_cs'] = [
        bool(leads.get((g, ct, c), set()) & ({canon(v, c) for v in m.split(',')} if c == 'pooled'
                                              else {v.upper() for v in m.split(',')}))
        for g, ct, c, m in zip(cs.gene, cs.cell_type, cs.cohort, cs.members)]

    for run in RUNS:
        for suffix in ('cs_id', 'jaccard', 'sv_shared', 'sv_led'):
            cs[f'{run}_{suffix}'] = None
    for _, g in cs.groupby(['gene', 'cell_type']):
        for r1, r2 in (('bioheart', 'tob'), ('pooled', 'bioheart'), ('pooled', 'tob')):
            x1, x2 = g[g.cohort == r1], g[g.cohort == r2]
            cand = sorted(((len(x.members_c & y.members_c) / len(x.members_c | y.members_c), i, j)
                           for i, x in x1.iterrows() for j, y in x2.iterrows() if x.members_c & y.members_c),
                          reverse=True)
            used = set()
            for jac, i, j in cand:
                if i in used or j in used:
                    continue
                used |= {i, j}
                shared = any(v.startswith('FED:') for v in cs.at[i, 'members_c'] & cs.at[j, 'members_c'])
                led = cs.at[i, 'lead_c'] == cs.at[j, 'lead_c'] and cs.at[i, 'lead_c'].startswith('FED:')
                for me, other, other_run in ((i, j, r2), (j, i, r1)):
                    cs.at[me, f'{other_run}_cs_id'] = cs.at[other, 'cs_id']
                    cs.at[me, f'{other_run}_jaccard'] = round(jac, 3)
                    cs.at[me, f'{other_run}_sv_shared'] = shared
                    cs.at[me, f'{other_run}_sv_led'] = led

    out = cs.drop(columns=['members_c', 'lead_c'])
    for run in RUNS:   # a run is never matched to itself
        out.loc[out.cohort == run, [f'{run}_{s}' for s in ('cs_id', 'jaccard', 'sv_shared', 'sv_led')]] = None
    out.to_csv(a.out, sep='\t', index=False)

    for run in RUNS:
        c = cs[cs.cohort == run]
        if c.empty:
            continue
        others = [r for r in RUNS if r != run and (cs.cohort == r).any()]
        msg = ', '.join(f'matched in {r} {c[f"{r}_cs_id"].notna().mean():.0%}' for r in others)
        print(f'  {run}: {len(c):,} CS ({msg}); with SV {(c.n_SV > 0).sum():,}, SV-led {(c.lead_type == "SV").sum():,}')
    print(f'wrote {a.out}')


if __name__ == '__main__':
    main()
