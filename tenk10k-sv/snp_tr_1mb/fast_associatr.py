#!/usr/bin/env python3
"""
associaTR-equivalent cis association, batched.

Keeps associaTR's own genotype loader (trtools.associaTR.load_and_filter_genotypes.load_trs), so locus
filters (no calls / one allele / non-major allele count < cutoff), allele lengths and per-locus called-sample
masks are identical by construction. Replaces only the per-locus statsmodels OLS with batched linear algebra:

    associaTR (6.0.1) fits, per locus, OLS(y_std ~ g_std + 1 + C_std) on called samples and reports
    coef/std(g)*std(y), se/std(g)*std(y), the OLS t-test p-value and the model R^2.

    Standardising y, g and C is an invertible affine transform, so the genotype coefficient, its SE and p are
    exactly those of OLS(y ~ g + 1 + C) on the raw scale. By Frisch-Waugh-Lovell that coefficient equals
    (g_r . y_r) / (g_r . g_r) with g_r, y_r the residuals after projecting out [1, C], and
    se^2 = RSS / (n - p) / (g_r . g_r), p = 2 * t_sf(|t|, n - p), p = number of regressors incl. g.

Loci whose called-sample set equals the full sample set share one projection (one QR for the whole window).
Loci with missing calls are solved individually on their own sample subset (exact, still numpy).

Output columns match associaTR's TSV: chrom, pos, alleles, n_samples_tested, locus_filtered, p_<name>,
coeff_<name>, se_<name>, regression_R^2, motif, period, ref_len, allele_frequency.

Usage (one gene, one or many pheno_cov files sharing the VCF window):
    fast_associatr.py --vcf X.vcf.bgz --region chr11:133412250-135412250 --vcftype eh \
        --pheno NK=NK.npy --pheno CD8_TEM=CD8_TEM.npy --out-prefix out/B3GAT1
"""
import argparse
import sys
import time

import cyvcf2
import numpy as np
from scipy import stats
from trtools.associaTR import load_and_filter_genotypes as lfg
import trtools.utils.tr_harmonizer as trh
from trtools.associaTR.associaTR import _merge_arrays

PVAL_PRECISION = 2
LAST_IDS = None       # VCF record IDs of the last run_window call, in output row order
EXACT_LOADER = False   # True: associaTR's load_trs per phenotype (slower; used to validate)
TIMING = dict(load=0.0, regress=0.0, n_full=0, n_partial=0)   # associaTR writes p with this many digits after the point in %e


def prepare_covars(all_samples, pheno_fname):
    """associaTR's covariate preparation (not --same-samples): returns sample filter, y, C (no intercept)."""
    covars = np.load(pheno_fname)
    covars = _merge_arrays(np.array(all_samples, dtype=float).reshape(-1, 1), covars)
    sample_filter = ~np.any(np.isnan(covars), axis=1)
    covars = covars[sample_filter, :]
    y = covars[:, 1].astype(float)
    C = covars[:, 2:].astype(float)
    return sample_filter, y, C


def ols_partial(g, called, p):
    """Exact OLS of y on [g, 1, C] restricted to the called samples, by downdating Z'Z for the few missing ones.

    Same algebra as the full-sample FWL path: with Z = [1, C] on the called set S,
      gg = g'g - (Z_S'g) (Z_S'Z_S)^-1 (Z_S'g),  gy = g'y_S - (Z_S'g) (Z_S'Z_S)^-1 (Z_S'y_S),
      yy = y_S'y_S - (Z_S'y_S) (Z_S'Z_S)^-1 (Z_S'y_S),
    where Z_S'Z_S = Z'Z - Z_M'Z_M and Z_S'y_S = Z'y - Z_M'y_M (M = missing samples).
    """
    miss = ~called
    Zm, ym = p['Z'][miss], p['y'][miss]
    ZtZ = p['ZtZ'] - Zm.T @ Zm
    Zty = p['Zty'] - Zm.T @ ym
    ys = p['y'][called]
    Ztg = p['Z'][called].T @ g
    L = np.linalg.solve(ZtZ, np.column_stack([Ztg, Zty]))
    gg = g @ g - Ztg @ L[:, 0]
    gy = g @ ys - Ztg @ L[:, 1]
    yy = ys @ ys - Zty @ L[:, 1]
    n_s = int(called.sum())
    b = gy / gg
    rss = yy - b * b * gg
    df = n_s - p['k'] - 1
    se = np.sqrt(rss / df / gg)
    tss = ys @ ys - ys.sum() ** 2 / n_s
    return b, se, 2 * stats.t.sf(abs(b / se), df), 1 - rss / tss


def load_window_once(vcf_fname, region, vcftype):
    """One pass over the window for ALL VCF samples, using associaTR's own harmonizer.

    Stores per locus the allele-index genotypes (int8, -1 = missing) and the unrounded allele lengths, so that the
    per-phenotype sample subsetting, non-major filter and allele-frequency string can be re-derived exactly as
    load_trs would for that phenotype's sample set.
    """
    vcf = cyvcf2.VCF(vcf_fname)
    vt = trh.InferVCFType(vcf, vcftype if vcftype else 'auto')
    region_start = int(region.split(':')[1].split('-')[0])
    loci = []
    for record in vcf(region):
        if record.POS < region_start:
            continue
        tr = trh.HarmonizeRecord(vcfrecord=record, vcftype=vt)
        lens = np.array([tr.ref_allele_length] + list(tr.alt_allele_lengths), dtype=float)
        idx = record.genotype.array()[:, :2].astype(np.int16)        # allele indices, -1 = missing
        loci.append(dict(chrom=tr.chrom, pos=tr.pos, id=record.ID, lens=lens, idx=idx, called=tr.GetCalledSamples(),
                         motif=tr.motif, ref_len=tr.ref_allele_length))
    return loci


def locus_for_samples(locus, samples, n_samples_total, non_major_cutoff):
    """Replicates load_trs for one locus and one boolean sample filter (associaTR not --beagle-dosages)."""
    called_samples_filter = locus['called'][samples]
    curr = samples & locus['called']
    n = int(curr.sum())
    len_alleles = [round(x, lfg.allele_len_precision) for x in locus['lens']]
    gts = locus['lens'][locus['idx'][curr]]                        # n x 2 lengths, unrounded (as GetLengthGenotypes)
    # GetAlleleFreqs(uselength=True) over curr samples, then clean_len_alleles
    if n:
        vals, counts = np.unique(gts.ravel(), return_counts=True)
        af = lfg.clean_len_alleles({v: c / (2 * n) for v, c in zip(vals, counts)})
    else:
        af = {}
    details = [locus['motif'], str(len(locus['motif'])), str(round(locus['ref_len'], lfg.allele_len_precision)),
               lfg.dict_str({key: '{:.2g}'.format(val) for key, val in af.items()})]
    if len(af) == 0:
        filt = 'No called samples'
    elif len(af) == 1:
        filt = 'Only one called allele'
    else:
        f = list(af.values()); f.pop(int(np.argmax(f)))
        filt = 'non-major allele count<{}'.format(non_major_cutoff) if np.sum(f) * n * 2 < non_major_cutoff else None
    return (None if filt else gts), np.unique(len_alleles), locus['chrom'], locus['pos'], called_samples_filter, filt, details


def run_window(vcf_fname, region, vcftype, phenos, non_major_cutoff=20):
    """phenos: dict name -> pheno_cov .npy path. Yields (name, row) per locus."""
    all_samples = cyvcf2.VCF(vcf_fname).samples
    prepared = {}
    for name, f in phenos.items():
        sf, y, C = prepare_covars(all_samples, f)
        # centre/scale covariates for conditioning only; the genotype coefficient is unchanged (affine reparam.)
        sd = C.std(0); sd[sd == 0] = 1
        Z = np.column_stack([np.ones(len(y)), (C - C.mean(0)) / sd])
        Q, _ = np.linalg.qr(Z)
        y_r = y - Q @ (Q.T @ y)
        prepared[name] = dict(sf=sf, y=y, Z=Z, ZtZ=Z.T @ Z, Zty=Z.T @ y, Q=Q, y_r=y_r,
                              tss=((y - y.mean()) ** 2).sum(), k=Z.shape[1])

    # all phenotypes for a gene must share a sample filter to share one genotype pass;
    # otherwise loop per phenotype (still one genotype pass each)
    groups = {}
    for name, p in prepared.items():
        groups.setdefault(p['sf'].tobytes(), []).append(name)

    out = {name: [] for name in phenos}
    t_load = time.time()
    loci = load_window_once(vcf_fname, region, vcftype) if not EXACT_LOADER else None
    global LAST_IDS
    LAST_IDS = [L['id'] for L in loci] if loci is not None else None   # VCF IDs, in output row order
    TIMING['load'] += time.time() - t_load
    for _, names in groups.items():
        sf = prepared[names[0]]['sf']
        if EXACT_LOADER:
            it = lfg.load_trs(vcf_fname, sf.copy(), region, non_major_cutoff, False, vcftype)
            next(it)   # field names
        else:
            it = (locus_for_samples(L, sf, len(sf), non_major_cutoff) for L in loci)
        full_G, full_meta, partial = [], [], []
        t_load = time.time()
        for idx, (gts, unique_alleles, chrom, pos, called, filt, details) in enumerate(it):
            alleles = ','.join(list(unique_alleles.astype(str)))
            n_called = int(np.sum(called))
            base = [idx, chrom, pos, alleles, n_called]   # idx = VCF record order, dropped on write
            if not filt and prepared[names[0]]['k'] + 1 >= n_called:
                filt = 'n covars >= n samples'
            if filt:
                for name in names:
                    out[name].append(base + [filt, np.nan, np.nan, np.nan, np.nan] + list(details))
                continue
            g = np.sum(gts, axis=1).astype(float)
            if called.all():
                full_G.append(g); full_meta.append((base, details))
            else:
                partial.append((g, called, base, details))

        TIMING['load'] += time.time() - t_load
        TIMING['n_full'] += len(full_G); TIMING['n_partial'] += len(partial)
        t_reg = time.time()
        for name in names:
            p = prepared[name]
            n, k = len(p['y']), p['k']
            if full_G:
                G = np.column_stack(full_G)
                G_r = G - p['Q'] @ (p['Q'].T @ G)
                gg = (G_r ** 2).sum(0)
                b = G_r.T @ p['y_r'] / gg
                rss = (p['y_r'] ** 2).sum() - b ** 2 * gg
                df = n - k - 1
                se = np.sqrt(rss / df / gg)
                pv = 2 * stats.t.sf(np.abs(b / se), df)
                r2 = 1 - rss / p['tss']
                for (base, details), bi, si, pi, ri in zip(full_meta, b, se, pv, r2):
                    out[name].append(base + ['False', pi, bi, si, ri] + list(details))
            for g, called, base, details in partial:
                bi, si, pi, ri = ols_partial(g, called, p)
                out[name].append(base + ['False', pi, bi, si, ri] + list(details))
        TIMING['regress'] += time.time() - t_reg
    for name in out:
        out[name].sort(key=lambda r: r[0])          # restore VCF record order, as associaTR writes it
        out[name] = [r[1:] for r in out[name]]
    return out


HEADER = ['chrom', 'pos', 'alleles', 'n_samples_tested', 'locus_filtered', 'p_{}', 'coeff_{}', 'se_{}',
          'regression_R^2', 'motif', 'period', 'ref_len', 'allele_frequency']


def write_tsv(path, name, rows):
    with open(path, 'w') as fh:
        fh.write('\t'.join(h.format(name) for h in HEADER) + '\n')
        for r in rows:
            r = list(r)
            if r[4] == 'False':
                r[5] = ('{:.' + str(PVAL_PRECISION) + 'e}').format(r[5])
            fh.write('\t'.join(str(x) for x in r) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--region', required=True)
    ap.add_argument('--vcftype', default='eh')
    ap.add_argument('--pheno', action='append', required=True, help='name=path/to/pheno_cov.npy (repeatable)')
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--non-major-cutoff', type=float, default=20)
    ap.add_argument('--exact-loader', action='store_true', help="use associaTR's load_trs per phenotype")
    a = ap.parse_args()
    global EXACT_LOADER
    EXACT_LOADER = a.exact_loader
    phenos = dict(x.split('=', 1) for x in a.pheno)
    t0 = time.time()
    res = run_window(a.vcf, a.region, a.vcftype, phenos, a.non_major_cutoff)
    for name, rows in res.items():
        write_tsv(f'{a.out_prefix}.{name}.tsv', name, rows)
    print(f'{len(next(iter(res.values())))} loci x {len(res)} phenotypes in {time.time() - t0:.1f}s '
          f'(genotype load {TIMING["load"]:.1f}s, regression {TIMING["regress"]:.1f}s; '
          f'{TIMING["n_full"]} complete-call loci, {TIMING["n_partial"]} with missing calls)', file=sys.stderr)


if __name__ == '__main__':
    main()
