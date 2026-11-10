#!/usr/bin/env python3
"""
Build per-cohort SuSiE inputs (SNP + TR + SV dosages, phenotype, covariates) for a batch of genes.

Runs inside the trtools image, next to a copy of ../snp_tr_1mb/fast_associatr.py, so genotype parsing and the
per-locus filters are associaTR's own (TRTools harmonizer, `load_trs` filters: no calls / one allele / non-major
allele count < 20), exactly as in the SNP/TR and SV eQTL scans.

Per gene and cohort it writes, under {out_dir}/{gene}/{cohort}/:

    X.f32.gz        float32 dosage matrix, samples x variants, row-major, NaN = missing call
    samples.txt     row order of X (numeric sample IDs, as in the VCFs)
    variants.tsv.gz one row per column of X: variant_id, variant_type (SNP/TR/SV), chrom, pos, end, motif,
                    svtype, overlapping_trs, and keep_<cell_type> (True if the locus passes associaTR's filters
                    in that cell type's sample set)
    pheno_<cell_type>.tsv.gz   sample, y, covariates (the SV eQTL pheno_cov array for that cohort)

and, when the manifest has `meta_set` (the FederateSV table), {out_dir}/{gene}/pooled_svs.tsv: the federated SV
pairs the pooled model may use (see write_pooled_svs).

Which variants are included -- the same windows and rules as the eQTL scans:
    SNP, TR: POS inside the gene's 1-Mb cis window (associaTR's rule, as in snp_tr_cis_runner.py)
    SV:      [POS, END] overlaps the window by >= 1 bp (associatr_runner_sv.py's bcftools pre-filter)

Two clean-ups the STR MCV pipeline also makes (str/fine-mapping/remove_STR_indels.py), done here from the VCF
records so they cover the whole 1-Mb window rather than only the 100-kb associaTR outputs:
    - indels whose inserted/deleted sequence is whole copies of a cyclic rotation of an overlapping TR's motif
      (or its reverse complement) are dropped -- they are the TR, called twice;
    - TRs duplicated at the same (chrom, pos, motif) keep their first record.
    One deliberate difference: every overlapping TR is checked, not only the first one found.

SVs that overlap a TR are KEPT and flagged in `overlapping_trs`. A VNTR expansion can be genotyped both as a TR and
as an SV INS/DEL; SuSiE will then split PIP between them, so a credible set holding an SV and an overlapping TR
should be read as one event.

Dosage = summed allele length per sample (associaTR's genotype), so a biallelic SV is its allele count and a
multiallelic CNV is copy number (sv_vcf_for_associatr.py encodes it as CN/0), both up to an affine shift that
residualisation removes.
"""
import gzip
import json
import os
import sys

import cyvcf2
import numpy as np
import pandas as pd
import trtools.utils.tr_harmonizer as trh

import fast_associatr as fa

NON_MAJOR_CUTOFF = 20
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(s):
    return s.translate(COMP)[::-1]


def rotations(s):
    return {s[i:] + s[:i] for i in range(len(s))}


def is_tr_indel(ref, alt, motif):
    """True if the indel's inserted/deleted bases are whole copies of a rotation of motif (or its revcomp)."""
    if len(ref) == len(alt) or not motif:
        return False
    seq = alt[1:] if len(alt) > len(ref) else ref[1:]
    if not seq or len(seq) % len(motif):
        return False
    rots = rotations(motif * (len(seq) // len(motif)))
    return seq in rots or revcomp(seq) in rots


def load_sv_chrom(vcf_fname):
    """Every SV record in a (per-chromosome) SV VCF, read end to end -- no index needed -- in
    fa.load_window_once's locus format, so fa.locus_for_samples applies associaTR's filters unchanged."""
    vcf = cyvcf2.VCF(vcf_fname)
    vt = trh.InferVCFType(vcf, 'eh')
    loci = []
    for record in vcf:
        tr = trh.HarmonizeRecord(vcfrecord=record, vcftype=vt)
        lens = np.array([tr.ref_allele_length] + list(tr.alt_allele_lengths), dtype=float)
        idx = record.genotype.array()[:, :2].astype(np.int16)
        loci.append(dict(chrom=tr.chrom, pos=tr.pos, end=int(record.INFO.get('END') or record.POS), id=record.ID,
                         lens=lens, idx=idx, called=tr.GetCalledSamples(), motif=tr.motif,
                         ref_len=tr.ref_allele_length, svtype=record.INFO.get('SVTYPE')))
    return loci


def svs_in_window(sv_loci, chrom, start, end):
    """associatr_runner_sv.py's rule: [POS, INFO/END] overlaps the window by >= 1 bp."""
    return [L for L in sv_loci if L['chrom'] == chrom and L['pos'] <= end and L['end'] >= start]


def dosage(locus, rows):
    """Summed allele length for the VCF sample indices `rows`; NaN where the call is missing."""
    idx = locus['idx'][rows]
    g = locus['lens'][np.clip(idx, 0, None)].sum(1)
    g[(idx < 0).any(1)] = np.nan
    return g


def tr_intervals(tr_loci):
    """TR spans (pos .. pos + ref_len * period) and motifs, for the indel and SV-overlap checks."""
    return [(L['pos'], L['pos'] + L['ref_len'] * len(L['motif']), L['motif']) for L in tr_loci]


def drop_tr_indels(snp_loci, trs):
    keep, n_drop = [], 0
    starts = np.array([t[0] for t in trs]) if trs else np.array([])
    for L in snp_loci:
        ref, alt = (L['motif'].split('-') + [''])[:2]    # mock-EH SNP VCF: RU = REF-ALT
        if len(ref) != len(alt) and len(starts):
            hits = [t for t in trs[:np.searchsorted(starts, L['pos'], side='right')] if t[1] >= L['pos']]
            if any(is_tr_indel(ref, alt, m) for _, _, m in hits):
                n_drop += 1
                continue
        keep.append(L)
    return keep, n_drop


def dedup_trs(tr_loci):
    seen, keep = set(), []
    for L in tr_loci:
        k = (L['chrom'], L['pos'], L['motif'])
        if k not in seen:
            seen.add(k)
            keep.append(L)
    return keep, len(tr_loci) - len(keep)


def load_pheno(path):
    """pheno_cov .npy: column 0 numeric sample ID, 1 phenotype, 2.. covariates. Rows with any NaN are dropped,
    as associaTR's covariate merge does."""
    a = np.load(path)
    a = a[~np.isnan(a).any(1)]
    return pd.DataFrame(a[:, 1:], index=a[:, 0].astype(np.int64).astype(str),
                        columns=['y'] + [f'c{i}' for i in range(1, a.shape[1] - 1)])


def build_gene(gene, man, snp_loci, tr_loci, snp_samples, tr_samples, sv_cache, out_dir, log):
    chrom, rng = gene['region'].split(':')
    start, end = map(int, rng.split('-'))
    trs = tr_intervals(tr_loci)
    tr_starts = np.array([t[0] for t in trs]) if trs else np.array([])
    for cohort, sv_vcf in man['sv_vcf'].items():
        phenos = {ct: load_pheno(by_coh[cohort]) for ct, by_coh in gene['pheno'].items() if cohort in by_coh}
        if not phenos:
            continue
        if cohort not in sv_cache:
            sv_cache[cohort] = (load_sv_chrom(sv_vcf), cyvcf2.VCF(sv_vcf).samples)
        sv_loci = svs_in_window(sv_cache[cohort][0], chrom, start, end)
        sv_samples = sv_cache[cohort][1]
        in_all = set(snp_samples) & set(tr_samples) & set(sv_samples)
        samples = sorted(set().union(*[set(p.index) for p in phenos.values()]) & in_all, key=int)
        n_absent = len(set().union(*[set(p.index) for p in phenos.values()])) - len(samples)
        pos_in = {name: {s: i for i, s in enumerate(lst)}
                  for name, lst in (('SNP', snp_samples), ('TR', tr_samples), ('SV', sv_samples))}

        cols, meta = [], []
        for vtype, loci in (('SNP', snp_loci), ('TR', tr_loci), ('SV', sv_loci)):
            rows = np.array([pos_in[vtype][s] for s in samples])
            n_vcf = len(pos_in[vtype])
            # per cell type: associaTR's locus filter on that cell type's sample set, in this VCF's sample order
            masks = {}
            for ct, p in phenos.items():
                m = np.zeros(n_vcf, bool)
                m[[pos_in[vtype][s] for s in p.index if s in pos_in[vtype] and s in in_all]] = True
                masks[ct] = m
            for L in loci:
                keep = {ct: fa.locus_for_samples(L, m, n_vcf, NON_MAJOR_CUTOFF)[5] is None for ct, m in masks.items()}
                if not any(keep.values()):
                    continue
                g = dosage(L, rows)
                if np.nanstd(g) == 0:
                    continue
                if vtype == 'SV':
                    hits = [f'{t[0]}:{t[2]}' for t in trs[:np.searchsorted(tr_starts, L['end'], side='right')]
                            if t[1] >= L['pos']] if trs else []
                    vid, motif, svtype, v_end, ovl = L['id'], '', L['svtype'], L['end'], ','.join(hits)
                elif vtype == 'TR':
                    vid, motif, svtype, v_end, ovl = f"{L['chrom']}:{L['pos']}:{L['motif']}", L['motif'], '', \
                        int(L['pos'] + L['ref_len'] * len(L['motif'])), ''
                else:
                    vid, motif, svtype, v_end, ovl = L['id'] or f"{L['chrom']}:{L['pos']}:{L['motif']}", \
                        L['motif'], '', L['pos'], ''
                cols.append(g.astype(np.float32))
                meta.append(dict(variant_id=vid, variant_type=vtype, chrom=L['chrom'], pos=L['pos'], end=v_end,
                                 motif=motif, svtype=svtype, overlapping_trs=ovl,
                                 **{f'keep_{ct}': k for ct, k in keep.items()}))

        d = os.path.join(out_dir, gene['gene'], cohort)
        os.makedirs(d, exist_ok=True)
        X = np.column_stack(cols) if cols else np.zeros((len(samples), 0), np.float32)
        with gzip.open(os.path.join(d, 'X.f32.gz'), 'wb', compresslevel=4) as fh:
            fh.write(np.ascontiguousarray(X, dtype='<f4').tobytes())
        with open(os.path.join(d, 'samples.txt'), 'w') as fh:
            fh.write('\n'.join(samples) + '\n')
        pd.DataFrame(meta).to_csv(os.path.join(d, 'variants.tsv.gz'), sep='\t', index=False)
        for ct, p in phenos.items():
            p.loc[[s for s in p.index if s in in_all]].rename_axis('sample').to_csv(
                os.path.join(d, f'pheno_{ct}.tsv.gz'), sep='\t')
        n_by = pd.Series([m['variant_type'] for m in meta]).value_counts().to_dict()
        print(f"{gene['gene']} {cohort}: {len(samples)} samples ({n_absent} pheno samples not in all VCFs), "
              f"variants {n_by}, {len(phenos)} cell types", flush=True)
        log.append(dict(gene=gene['gene'], cohort=cohort, n_samples=len(samples), n_pheno_not_in_vcfs=n_absent,
                        **{f'n_{k}': v for k, v in n_by.items()}))


POOLABLE = lambda c: c.startswith('concordant_') or c == 'dup_cnv'   # meta_runner_sv.py's pooled classes


def load_pairs(meta_set_path):
    """FederateSV pairing table, filtered exactly as meta_runner_sv._load_pairs (include_dup_cnv=True): keep
    concordant_* and dup_cnv (del_cnv never: DEL dosage runs opposite to CN), first VARID of any comma-joined
    cell, keys upper-cased because associaTR upper-cases the VARID it embeds while the VCF keeps mixed case.

    Pass the table the meta-analysis actually pooled:
    - 7pc/v2_F_filtered: meta_analysis_set_v2_fixed_F_filtered.csv (15,109 pairs; holds all 14,424 tested pairs);
    - 7pc/v1_F_filtered: the original meta_analysis_set.csv restricted to pairs passing the F filter
      (meta_analysis_set.as_used_v1_F_filtered.csv): its 12,972 tested pairs all lie in that intersection, and none
      of the pairs added by meta_analysis_set_fixed.csv were tested.

    A pair's identity is its (BioHEART VARID, TOB VARID); federate_id is NOT unique in the fixed tables, so the
    pooled column ID (pair_id) is the federate_id, suffixed with the BioHEART VARID wherever it repeats.
    """
    pairs = pd.read_csv(meta_set_path)
    pairs = pairs[pairs['pair_class'].map(POOLABLE)].copy()
    for col in ('bioheart_vids', 'tob_vids'):
        pairs[col] = pairs[col].astype(str).str.split(',').str[0].str.strip()
    pairs['bh_key'] = pairs['bioheart_vids'].str.upper()
    pairs['tob_key'] = pairs['tob_vids'].str.upper()
    for key in ('bh_key', 'tob_key'):
        if pairs[key].duplicated().any():
            raise ValueError(f'{meta_set_path}: a {key} VARID appears in more than one pair; pairing is not 1:1')
    rep = pairs['federate_id'].duplicated(keep=False)
    pairs['pair_id'] = pairs['federate_id'].where(~rep, pairs['federate_id'] + '#' + pairs['bioheart_vids'])
    return pairs


def write_pooled_svs(gene, out_dir, pairs):
    """{gene}/pooled_svs.tsv: one row per federated pair with both SVs in their cohort's window table.

    Whether each side passes associaTR's filters in a given cell type is decided later, per cell type, from the
    cohort's keep_<cell_type> column (the meta-analysis needs the SV tested in both cohorts). SV IDs that occur
    more than once in a cohort's window are dropped, as meta_runner_sv.py's drop_duplicates(keep=False) does.
    """
    d = os.path.join(out_dir, gene)
    tabs = {}
    for coh in ('bioheart', 'tob'):
        f = os.path.join(d, coh, 'variants.tsv.gz')
        if not os.path.exists(f):
            return 0
        v = pd.read_csv(f, sep='\t', usecols=['variant_id', 'variant_type'])
        sv = v.loc[v.variant_type == 'SV', 'variant_id'].str.upper()
        tabs[coh] = set(sv[~sv.duplicated(keep=False)])
    m = pairs[pairs.bh_key.isin(tabs['bioheart']) & pairs.tob_key.isin(tabs['tob'])]
    m[['pair_id', 'federate_id', 'pair_class', 'bioheart_vids', 'tob_vids', 'bioheart_svtype', 'tob_svtype']].rename(
        columns={'bioheart_vids': 'bioheart_id', 'tob_vids': 'tob_id'}).to_csv(
        os.path.join(d, 'pooled_svs.tsv'), sep='\t', index=False)
    return len(m)


def main(manifest_path, out_dir):
    man = json.loads(os.path.expandvars(open(manifest_path).read()))
    snp_samples = cyvcf2.VCF(man['snp_vcf']).samples
    tr_samples = cyvcf2.VCF(man['tr_vcf']).samples
    log, sv_cache = [], {}
    pairs = load_pairs(man['meta_set']) if man.get('meta_set') else None
    for gene in man['genes']:
        snp_loci = fa.load_window_once(man['snp_vcf'], gene['region'], 'eh')
        tr_loci, n_dup = dedup_trs(fa.load_window_once(man['tr_vcf'], gene['region'], 'eh'))
        snp_loci, n_tr_indel = drop_tr_indels(snp_loci, tr_intervals(tr_loci))
        print(f"{gene['gene']}: dropped {n_tr_indel} TR-like indels, {n_dup} duplicate TRs", flush=True)
        build_gene(gene, man, snp_loci, tr_loci, snp_samples, tr_samples, sv_cache, out_dir, log)
        if pairs is not None:
            n = write_pooled_svs(gene['gene'], out_dir, pairs)
            print(f"{gene['gene']}: {n} federated SV pairs with both sides in the window", flush=True)
    pd.DataFrame(log).to_csv(os.path.join(out_dir, 'prep_log.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
