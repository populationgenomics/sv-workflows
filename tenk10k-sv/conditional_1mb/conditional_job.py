#!/usr/bin/env python3
# pylint: disable=too-many-locals,too-many-statements
"""
Conditional cis association for SV-led eGene x cell type pairs. Runs inside the trtools image, next to copies of
../snp_tr_1mb/fast_associatr.py, ../snp_tr_1mb/meta_np.py and ../finemap_mcv/finemap_prep.py.

For each pair, every SNP, TR and SV in the gene's 1-Mb cis window is tested three times, per cohort, then
meta-analysed:
    marginal          no extra covariate (the "before"; same code, samples and covariates as the scans)
    cond_lead_snptr   + dosage of the lead SNP/TR (largest |z| among SNPs and TRs in the 1-Mb re-test)
    cond_lead_sv      + dosage of the lead SV (the cohort's own record of the federated pair)

Conditioning = the conditioning variant's dosage (associaTR's summed allele length, the regressor associaTR itself
used) is appended as one more column of the cohort's pheno_cov array, so fast_associatr.run_window fits
OLS(y ~ g + 1 + C + g_cond) with associaTR's own locus filters. Missing conditioning calls are mean-imputed within
the cohort so the sample set is identical before and after conditioning (n imputed is reported).

Variant windows, exactly as in the eQTL scans: SNPs/TRs by POS (fast_associatr.load_window_once, i.e. the 1-Mb
re-test); SVs by [POS, END] overlap (finemap_prep.svs_in_window, i.e. associatr_runner_sv.py), read from the
per-cohort SV VCF end to end. Meta-analysis: SNPs/TRs by record (jointly called, same VCF); SVs by federated pair
from the pairing table (concordant_* and dup_cnv), as meta_runner_sv.py.

Usage: python3 conditional_job.py manifest.json assoc.tsv.gz summary.tsv
"""
import gzip
import json
import os
import sys

import cyvcf2
import numpy as np
import pandas as pd

import fast_associatr as fa
import finemap_prep as fp
from meta_np import metagen_dl

COLS = ['chrom', 'pos', 'alleles', 'n_samples_tested', 'locus_filtered', 'p', 'coeff', 'se', 'r2',
        'motif', 'period', 'ref_len', 'allele_frequency']
COHORT_SFX = {'bioheart': '1', 'tob': '2'}
CONDITIONS = ('marginal', 'cond_lead_snptr', 'cond_lead_sv')
GWS = 5e-8

# ---- serve preloaded windows to fast_associatr.run_window (one genotype pass per gene, reused for every run) ----
WINDOWS = {}
_ORIG_LOADER = fa.load_window_once


def _cached_loader(vcf, region, vcftype):
    key = (vcf, region)
    if key not in WINDOWS:
        WINDOWS[key] = _ORIG_LOADER(vcf, region, vcftype)
    return WINDOWS[key]


fa.load_window_once = _cached_loader
SAMPLES = {}


def vcf_samples(vcf):
    if vcf not in SAMPLES:
        SAMPLES[vcf] = cyvcf2.VCF(vcf).samples
    return SAMPLES[vcf]


def sample_rows(vcf, sample_ids):
    """VCF column index for each pheno sample id (-1 if absent)."""
    pos = {float(s): i for i, s in enumerate(vcf_samples(vcf))}
    return np.array([pos.get(float(s), -1) for s in sample_ids])


def dosage_on(locus, vcf, sample_ids):
    """Summed allele length per pheno sample; NaN if missing or absent from the VCF."""
    rows = sample_rows(vcf, sample_ids)
    g = np.full(len(rows), np.nan)
    ok = rows >= 0
    if ok.any():
        g[ok] = fp.dosage(locus, rows[ok])
    return g


def load_pheno(path):
    a = np.load(path, allow_pickle=True).astype(float)
    return a


def write_cond_pheno(pheno, g_cond, path):
    """pheno_cov with the conditioning dosage appended as the last covariate (mean-imputed within cohort)."""
    g = g_cond.copy()
    miss = np.isnan(g)
    if miss.all():
        raise ValueError('conditioning variant has no called samples')
    g[miss] = np.nanmean(g)
    np.save(path, np.column_stack([pheno, g]))
    return int(miss.sum())


def r2_with(G, x):
    """Squared Pearson r of each column of G with x, mean-imputing missing values (G: samples x loci)."""
    G = np.where(np.isnan(G), np.nanmean(G, axis=0, keepdims=True), G)
    x = np.where(np.isnan(x), np.nanmean(x), x)
    Gc, xc = G - G.mean(0), x - x.mean()
    den = np.sqrt((Gc ** 2).sum(0) * (xc ** 2).sum())
    with np.errstate(invalid='ignore', divide='ignore'):
        return ((Gc.T @ xc) / den) ** 2


def run(vcf, region, phenos):
    res = fa.run_window(vcf, region, 'eh', phenos)
    ids = list(fa.LAST_IDS)
    out = {}
    for name, rows in res.items():
        d = pd.DataFrame(rows, columns=COLS)
        d.insert(0, 'id', ids)
        d['ok'] = (d.locus_filtered.astype(str) == 'False') & np.isfinite(d.coeff.astype(float)) & \
            np.isfinite(d.se.astype(float)) & (d.se.astype(float) > 0)
        out[name] = d
    return out


def meta_rows(a, b, keys):
    """a, b: per-cohort frames aligned row by row; keys: extra columns for the output."""
    ok = a.ok.to_numpy() & b.ok.to_numpy()
    if not ok.any():
        return pd.DataFrame()
    a, b = a[ok].reset_index(drop=True), b[ok].reset_index(drop=True)
    m = metagen_dl(np.column_stack([a.coeff.astype(float), b.coeff.astype(float)]),
                   np.column_stack([a.se.astype(float), b.se.astype(float)]))
    cols = {k: (v[ok].reset_index(drop=True) if isinstance(v, pd.Series) else v) for k, v in keys.items()}
    cols.update({
        'n_1': a.n_samples_tested, 'n_2': b.n_samples_tested,
        'coeff_1': a.coeff.astype(float), 'se_1': a.se.astype(float), 'pval_1': a.p.astype(float),
        'coeff_2': b.coeff.astype(float), 'se_2': b.se.astype(float), 'pval_2': b.p.astype(float),
        'coeff_meta_fixed': m['coeff_meta_fixed'], 'se_meta_fixed': m['se_meta_fixed'],
        'pval_meta_fixed': m['pval_meta_fixed'], 'i2_meta': m['i2_meta'], 'pval_q_meta': m['pval_q_meta']})
    return pd.DataFrame(cols)


def main(manifest_path, assoc_path, summary_path):
    man = json.loads(os.path.expandvars(open(manifest_path).read()))   # ${BATCH_TMPDIR} paths from a quoted heredoc
    pairs_tab = fp.load_pairs(man['meta_set'])
    tob_of = dict(zip(pairs_tab.bh_key, pairs_tab.tob_key))
    class_of = dict(zip(pairs_tab.bh_key, pairs_tab.pair_class))
    sv_chrom = {}
    tmp = os.path.abspath('cond_pheno')
    os.makedirs(tmp, exist_ok=True)
    assoc_out = gzip.open(assoc_path, 'wt')
    header = True
    summary = []

    for g in man['genes']:
        gene, chrom, region = g['gene'], g['chrom'], g['region']
        start, end = (int(x) for x in region.split(':')[1].split('-'))
        WINDOWS.clear()
        snp_loci = _cached_loader(man['snp_vcf'], region, 'eh')
        tr_loci = _cached_loader(man['tr_vcf'], region, 'eh')
        sv_win = {}
        for coh in COHORT_SFX:
            vcf = man['sv_vcf'][coh]
            if vcf not in sv_chrom:
                sv_chrom.clear()   # one chromosome per job, but stay bounded
                sv_chrom[vcf] = fp.load_sv_chrom(vcf)
            sv_win[coh] = fp.svs_in_window(sv_chrom[vcf], chrom, start, end)
            WINDOWS[(vcf, region)] = sv_win[coh]
        nonsv = {('SNP', man['snp_vcf']): snp_loci, ('TR', man['tr_vcf']): tr_loci}

        for pr in g['pairs']:
            ct = pr['cell_type']
            note = []
            ph = {coh: load_pheno(pr['pheno'][coh]) for coh in COHORT_SFX}
            # ---- conditioning variants ----
            lt, lid = pr['lead_nonsv']['variant_type'], pr['lead_nonsv']['id']
            lvcf = man['snp_vcf'] if lt == 'SNP' else man['tr_vcf']
            lead_loc = [L for L in nonsv[(lt, lvcf)] if str(L['id']) == str(lid)]
            sv_ids = {'bioheart': pr['lead_sv']['bioheart_id'], 'tob': pr['lead_sv']['tob_id']}
            sv_loc = {coh: [L for L in sv_win[coh] if str(L['id']).upper() == str(sv_ids[coh]).upper()]
                      for coh in COHORT_SFX}
            if len(lead_loc) != 1:
                note.append(f'lead {lt} {lid}: {len(lead_loc)} records in window')
            for coh in COHORT_SFX:
                if len(sv_loc[coh]) != 1:
                    note.append(f'lead SV {coh} {sv_ids[coh]}: {len(sv_loc[coh])} records in window')
            if note:
                summary.append(dict(gene_name=gene, cell_type=ct, chunk=man.get('chunk', ''), note='; '.join(note)))
                print(f'{gene} {ct}: SKIPPED ({"; ".join(note)})', flush=True)
                continue
            g_lead = {coh: dosage_on(lead_loc[0], lvcf, ph[coh][:, 0]) for coh in COHORT_SFX}
            g_sv = {coh: dosage_on(sv_loc[coh][0], man['sv_vcf'][coh], ph[coh][:, 0]) for coh in COHORT_SFX}

            # ---- LD of every variant with the lead SV and the lead SNP/TR, per cohort, on this cell type's samples ----
            ld = {}
            for coh, sfx in COHORT_SFX.items():
                for (vt, vcf), loci in list(nonsv.items()) + [(('SV', man['sv_vcf'][coh]), sv_win[coh])]:
                    if not loci:
                        continue
                    rows = sample_rows(vcf, ph[coh][:, 0])
                    okr = rows >= 0
                    G = np.full((len(rows), len(loci)), np.nan)
                    for j, L in enumerate(loci):
                        G[okr, j] = fp.dosage(L, rows[okr])
                    for tag, x in (('sv', g_sv[coh]), ('lead', g_lead[coh])):
                        for L, r2 in zip(loci, r2_with(G, x)):
                            ld[(vt, str(L['id']).upper(), sfx, tag)] = r2

            n_imp = {}
            for cond in CONDITIONS:
                paths = {}
                for coh in COHORT_SFX:
                    if cond == 'marginal':
                        paths[coh] = pr['pheno'][coh]
                    else:
                        x = g_lead[coh] if cond == 'cond_lead_snptr' else g_sv[coh]
                        paths[coh] = os.path.join(tmp, f'{gene}_{ct}_{cond}_{coh}.npy')
                        n_imp[(cond, coh)] = write_cond_pheno(ph[coh], x, paths[coh])
                cond_on = {'marginal': '', 'cond_lead_snptr': f'{lt}:{lid}',
                           'cond_lead_sv': f"SV:{sv_ids['bioheart']}|{sv_ids['tob']}"}[cond]
                parts = []
                # SNPs and TRs: both cohorts in one VCF, rows aligned
                for (vt, vcf) in nonsv:
                    r = run(vcf, region, {coh: paths[coh] for coh in COHORT_SFX})
                    a, b = r['bioheart'], r['tob']
                    if cond == 'cond_lead_snptr' and vt == lt:   # the conditioning variant itself is collinear
                        self_ = a.id.astype(str) == str(lid)
                        a.loc[self_, 'ok'] = False; b.loc[self_, 'ok'] = False
                    parts.append(meta_rows(a, b, {
                        'variant_type': vt, 'chr': a.chrom, 'pos': a.pos.astype(int), 'end': a.pos.astype(int),
                        'id': a.id.astype(str), 'tob_id': '', 'pair_class': '', 'motif': a.motif}))
                # SVs: one VCF per cohort, paired through the federated table
                sv = {}
                for coh in COHORT_SFX:
                    vcf = man['sv_vcf'][coh]
                    d = run(vcf, region, {coh: paths[coh]})[coh]
                    d['end'] = [L['end'] for L in sv_win[coh]]
                    d['key'] = d.id.astype(str).str.upper()
                    if cond == 'cond_lead_sv':
                        d.loc[d.key == str(sv_ids[coh]).upper(), 'ok'] = False
                    sv[coh] = d.drop_duplicates('key', keep=False)   # as meta_runner_sv.py
                a = sv['bioheart'].assign(tob_key=sv['bioheart'].key.map(tob_of))
                a = a[a.tob_key.notna()]
                b = sv['tob'].set_index('key').reindex(a.tob_key).reset_index()
                hit = b.id.notna().to_numpy()
                a, b = a[hit].reset_index(drop=True), b[hit].reset_index(drop=True)
                if len(a):
                    parts.append(meta_rows(a, b, {
                        'variant_type': 'SV', 'chr': a.chrom, 'pos': a.pos.astype(int), 'end': a.end.astype(int),
                        'id': a.id.astype(str), 'tob_id': b.id.astype(str), 'pair_class': a.key.map(class_of),
                        'motif': ''}))
                res = pd.concat([p for p in parts if len(p)], ignore_index=True)
                for sfx in COHORT_SFX.values():
                    # SVs are keyed by each cohort's own record ID; SNPs/TRs share one ID
                    ids = res.id if sfx == '1' else res.tob_id.where(res.variant_type == 'SV', res.id)
                    for tag in ('sv', 'lead'):
                        res[f'r2_{tag}_{sfx}'] = [ld.get((vt, str(i).upper(), sfx, tag), np.nan)
                                                  for vt, i in zip(res.variant_type, ids)]
                res.insert(0, 'cond_on', cond_on)
                res.insert(0, 'condition', cond)
                res.insert(0, 'cell_type', ct)
                res.insert(0, 'gene_name', gene)
                res.to_csv(assoc_out, sep='\t', index=False, header=header, float_format='%.6g')
                header = False
                pr.setdefault('_res', {})[cond] = res

            # ---- one summary row per pair ----
            R = pr.pop('_res')
            z = lambda d: d.coeff_meta_fixed / d.se_meta_fixed

            def stat(cond, vt, vid):
                d = R[cond]
                d = d[(d.variant_type == vt) & (d.id.astype(str).str.upper() == str(vid).upper())]
                if d.empty:
                    return np.nan, np.nan
                return float(z(d).iloc[0]), float(d.pval_meta_fixed.iloc[0])

            def top(cond):
                d = R[cond].assign(absz=lambda x: z(x).abs()).sort_values('absz', ascending=False)
                r = d.iloc[0]
                return f"{r.variant_type}:{r.id}", float(r.absz), int((d.pval_meta_fixed < GWS).sum())

            row = dict(gene_name=gene, cell_type=ct, chrom=chrom, region=region, chunk=man.get('chunk', ''),
                       lead_sv_bioheart=sv_ids['bioheart'], lead_sv_tob=sv_ids['tob'],
                       lead_nonsv_type=lt, lead_nonsv_id=lid, input_absz_sv=pr['lead_sv']['absz'],
                       input_absz_nonsv=pr['lead_nonsv']['absz'], note='')
            for cond in CONDITIONS:
                row[f'z_sv_{cond}'], row[f'p_sv_{cond}'] = stat(cond, 'SV', sv_ids['bioheart'])
                row[f'z_lead_{cond}'], row[f'p_lead_{cond}'] = stat(cond, lt, lid)
                row[f'top_{cond}'], row[f'top_absz_{cond}'], row[f'n_gws_{cond}'] = top(cond)
            for coh, sfx in COHORT_SFX.items():
                row[f'r2_sv_lead_{sfx}'] = ld.get((lt, str(lid).upper(), sfx, 'sv'), np.nan)
                row[f'n_imputed_lead_{sfx}'] = n_imp.get(('cond_lead_snptr', coh))
                row[f'n_imputed_sv_{sfx}'] = n_imp.get(('cond_lead_sv', coh))
            summary.append(row)
            print(f"{gene} {ct}: |z| SV {abs(row['z_sv_marginal']):.1f} -> {abs(row['z_sv_cond_lead_snptr']):.1f} "
                  f"given {lt}; {lt} {abs(row['z_lead_marginal']):.1f} -> {abs(row['z_lead_cond_lead_sv']):.1f} "
                  f"given SV", flush=True)
    assoc_out.close()
    pd.DataFrame(summary).to_csv(summary_path, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:4])
