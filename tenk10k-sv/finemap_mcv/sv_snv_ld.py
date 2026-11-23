#!/usr/bin/env python3
"""
Best-tagging SNV for each lead sc-eSV on one chromosome: r^2 between the SV and every small variant within +-1 Mb,
computed the way the pooled SuSiE model sees the data.

    python3 sv_snv_ld.py <snp_vcf> <bioheart_sv_vcf> <tob_sv_vcf> <leads.tsv> <out.tsv> [window_bp=1000000]

<leads.tsv>: one row per lead SV pair (columns bioheart_vid, tob_vid), all on this chromosome.

Genotypes
    SNV / indel  the jointly called common-variant VCF (BioHEART + TOB samples). Its ALT is a symbolic placeholder
                 (<STR1>), so the alleles are read from the ID (chr:pos:REF:ALT); a variant is an SNV when both
                 alleles are one base. Dosage = alt-allele count from GT.
    SV           each cohort's associaTR SV VCF; dosage = sum of the GT allele indices (alt count for biallelic SVs,
                 copy number for CNVs), as in the eQTL scan.
    Pooled column: the BioHEART SV's dosage for BioHEART samples and the TOB SV's for TOB samples, as in the pooled
    fine-mapping. Missing calls are mean-imputed and every column is centred WITHIN each cohort before stacking, so
    cohort mean differences (e.g. DUP vs CNV encoding in dup_cnv pairs) do not create LD.

Output, one row per lead SV pair:
    status            ok / SV not found / no variants or constant SV
    n_variants        small variants in the window
    best_r2, best_variant, best_variant_pos              best tag among all small variants (SNVs and indels)
    best_snv_r2, best_snv                                best tag among SNVs only
    n_r2_ge_095, n_r2_ge_08, n_snv_r2_ge_095             counts of tags above thresholds
    sv_pos, sv_end                                       BioHEART SV interval used for the window
"""
import sys

import cyvcf2
import numpy as np
import pandas as pd


def sv_dosages(path, ids):
    v = cyvcf2.VCF(path)
    out = {}
    for r in v:
        if r.ID in ids:
            g = np.array([np.nan if any(a < 0 for a in x[:2]) else x[0] + x[1] for x in r.genotypes], dtype=float)
            out[r.ID] = (g, r.POS, max(r.INFO.get('END') or r.POS, r.POS))
    return v.samples, out


def centre(m, n_first):
    """Mean-impute and centre each row within the first n_first columns (BioHEART) and the rest (TOB)."""
    m = m.astype(np.float32)
    m[m < 0] = np.nan
    for sl in (slice(0, n_first), slice(n_first, None)):
        with np.errstate(invalid='ignore'):
            mu = np.nanmean(m[:, sl], axis=1, keepdims=True)
        m[:, sl] = np.nan_to_num(m[:, sl] - np.nan_to_num(mu))
    return m


def main(snp_vcf, bh_vcf, tob_vcf, leads_tsv, out_tsv, window=1_000_000):
    leads = pd.read_csv(leads_tsv, sep='\t').drop_duplicates(['bioheart_vid', 'tob_vid'])
    bh_s, bh = sv_dosages(bh_vcf, set(leads.bioheart_vid))
    tb_s, tb = sv_dosages(tob_vcf, set(leads.tob_vid))

    v = cyvcf2.VCF(snp_vcf)
    idx = {s: i for i, s in enumerate(v.samples)}
    bh_keep = [k for k, s in enumerate(bh_s) if s in idx]
    tb_keep = [k for k, s in enumerate(tb_s) if s in idx]
    cols = np.array([idx[bh_s[k]] for k in bh_keep] + [idx[tb_s[k]] for k in tb_keep])
    nb = len(bh_keep)
    print(f'samples: SNV VCF {len(v.samples)} | BioHEART SV {len(bh_s)} matched {nb} | TOB SV {len(tb_s)} '
          f'matched {len(tb_keep)} | lead SV pairs {len(leads)}', flush=True)

    pos, ids, snv, rows_g = [], [], [], []
    for r in v:
        gt = r.gt_types                                  # 0 hom-ref, 1 het, 3 hom-alt, 2 unknown
        d = np.where(gt == 3, 2, gt).astype(np.int8)
        d[gt == 2] = -1
        vid = r.ID or f'{r.CHROM}:{r.POS}:{r.REF}:{r.ALT[0] if r.ALT else ""}'
        al = vid.split(':')[2:]
        pos.append(r.POS)
        ids.append(vid)
        snv.append(len(al) == 2 and len(al[0]) == 1 and len(al[1]) == 1)
        rows_g.append(d[cols])
    pos, snv, G = np.array(pos), np.array(snv), np.vstack(rows_g)
    del rows_g
    print(f'small variants {len(pos):,}, of which SNVs {int(snv.sum()):,}', flush=True)

    rows = []
    for r in leads.itertuples():
        base = dict(bioheart_vid=r.bioheart_vid, tob_vid=r.tob_vid)
        if r.bioheart_vid not in bh or r.tob_vid not in tb:
            rows.append({**base, 'status': 'SV not found'})
            continue
        x = np.concatenate([bh[r.bioheart_vid][0][bh_keep], tb[r.tob_vid][0][tb_keep]])
        p0, p1 = bh[r.bioheart_vid][1], bh[r.bioheart_vid][2]
        w = np.where((pos >= p0 - window) & (pos <= p1 + window))[0]
        xc = centre(x[None, :], nb)[0]
        if len(w) == 0 or float(xc @ xc) == 0:
            rows.append({**base, 'status': 'no variants or constant SV'})
            continue
        s = centre(G[w], nb)
        den = np.sqrt((s * s).sum(1) * float(xc @ xc))
        den[den == 0] = np.inf
        r2 = ((s @ xc) / den) ** 2
        k = int(np.argmax(r2))
        m = snv[w]
        ks = int(np.argmax(np.where(m, r2, -1))) if m.any() else None
        rows.append({**base, 'status': 'ok', 'n_variants': len(w),
                     'best_r2': float(r2[k]), 'best_variant': ids[w[k]], 'best_variant_pos': int(pos[w[k]]),
                     'best_snv_r2': float(r2[ks]) if ks is not None else np.nan,
                     'best_snv': ids[w[ks]] if ks is not None else None,
                     'n_r2_ge_095': int((r2 >= 0.95).sum()), 'n_r2_ge_08': int((r2 >= 0.8).sum()),
                     'n_snv_r2_ge_095': int(((r2 >= 0.95) & m).sum()), 'sv_pos': p0, 'sv_end': p1})
    pd.DataFrame(rows).to_csv(out_tsv, sep='\t', index=False)
    print(f'wrote {len(rows)} rows', flush=True)


if __name__ == '__main__':
    a = sys.argv[1:]
    main(*a[:5], window=int(a[5]) if len(a) > 5 else 1_000_000)
