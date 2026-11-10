#!/usr/bin/env python3
# pylint: disable=too-many-locals,too-many-statements,too-many-branches
"""
Cell-type specificity of eSVs: the five pairwise scenarios of the TenK10K flagship paper (Cuomo et al. 2025,
Methods "cell-type specificity"), as also used for the TR paper (str/cell-type-spec/).

Unit: every eQTL = the lead eSV of an eGene in its "primary" cell type A (gene-level FDR < 5%), keeping exact |z| ties.
For each of the other 27 cell types B, one pairwise scenario:

    1  eGene expressed only in the primary cell type          (gene not on B's expressed-gene list)
    2  expressed in both, eQTL only in the primary cell type  (expressed in B, gene-level q >= 0.05 in B)
    3  independent eQTL variants for the eGene in each        (eGene in B; B's lead SV(s): max r2 with A's lead < 0.5)
    4  same eQTL variant, or variants in LD (r2 >= 0.5), with concordant direction of effect
    5  same eQTL variant, or variants in LD (r2 >= 0.5), with discordant direction of effect

Rescue of scenario 2 (power in rarer cell types): A's lead SV is looked up in B, and a two-study random-effects
meta-analysis (REML, as metafor::rma, the flagship's v4.0.0 default) of its effect in A and in B is compared with
A's p-value. A >= 1e5-fold reduction (p_meta / p_A < 1e-5) reclassifies B as 4 or 5 by the concordance of the effect
direction. Compared on the log scale, so it works when p underflows.
    --rescue-effect abs (default): meta-analysis of |beta_A| and |beta_B|. The flagship text reclassifies rescued pairs
        as 4 OR 5 by concordance, which needs pooled magnitudes (with signed effects, discordant pairs can never lower
        p); the TR implementation of the same procedure (meta_eqtls_runner.py) pools abs(coeff).
    --rescue-effect signed: signed effects (sensitivity; rescues can then only be scenario 4).

Direction of effect for two different SVs is aligned through the sign of their genotype correlation r: effects are
concordant when sign(beta_A) == sign(beta_B) * sign(r). (A positive beta for an SV whose alt allele is anti-correlated
with A's alt allele is the same allelic direction as a negative beta for A.)

Expression (scenario 1 vs 2): B's expressed-gene list (scRNA_gene_lists/1_min_pct_cells_expressed), in either cohort.
SV identity: the federated pair (BioHEART VARID | TOB VARID). LD between two SVs: Pearson r of associaTR's dosage
(summed allele length) over all samples of each cohort's SV VCF, within cohort, combined across cohorts weighted by the
number of samples with both genotypes (within-cohort r is unaffected by the DUP/CNV coding offset of dup_cnv pairs).
Per eGene, the aggregate score is the maximum scenario over all pairs of cell types (flagship Fig. 4a).

Needs trtools (SV dosages via ../finemap_mcv/finemap_prep.py) and pyarrow.

Outputs in --out-dir:
    esv_cell_type_specificity.tsv   one row per eQTL (lead eSV x primary cell type): scenario_<cell type> (28 columns),
                                    rescued_cts, specific_per_row, n_ct_shared, max_scenario (per eGene), cell_spec_num
    sv_pair_ld.tsv                  SV-SV r and r2 used (per cohort and combined)
    rescue_tests.tsv                every scenario-2 meta-analysis: beta/se in A and B, tau2, p_meta, rescued yes/no

python3 sv_cell_type_specificity.py --sumstats full_sumstats.parquet --out-dir out/          (local; reads GCS with gcloud)
On Hail Batch: cell_type_spec_runner.py (gather job + this script in the trtools image on local inputs).
"""
import io
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'finemap_mcv'))
sys.path.insert(0, os.path.join(HERE, '..', 'snp_tr_1mb'))

CELL_TYPES = ('CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,'
              'CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,'
              'pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC').split(',')
V2 = 'gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v2_F_filtered'
INPUTS = 'gs://cpg-tenk10k-test/sv/input_files'
GENE_LISTS = INPUTS + '/{cohort}/7pc/scRNA_gene_lists/1_min_pct_cells_expressed/{ct}'
SV_VCF = ('gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/'
          '{cohort}_common_chr_specific/hail_filtered_{chrom}.vcf.bgz')
COHORT_DIRS = {'bioheart': 'bioheart_n968', 'tob': 'tob_n935'}
COHORTS = tuple(COHORT_DIRS)
LD_R2 = 0.5            # flagship: independent if r2 < 0.5
RESCUE_LOG10 = -5      # flagship: p_meta / p_primary < 1e-5
TIE = 1e-6             # |z| ties for the lead


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, executable='/bin/bash').stdout


def log10p(z):
    """log10 of the two-sided normal p, without underflow."""
    return (stats.norm.logsf(np.abs(z)) + np.log(2)) / np.log(10)


def rma_reml(y, v, max_iter=100, tol=1e-10):
    """Random-effects meta-analysis, REML tau2 by Fisher scoring (metafor::rma defaults: method='REML', test='z').
    Returns (estimate, se, tau2, z)."""
    y, v = np.asarray(y, float), np.asarray(v, float)
    tau2 = max(0.0, float(np.var(y, ddof=1) - v.mean())) if len(y) > 1 else 0.0
    for _ in range(max_iter):
        w = 1 / (v + tau2)
        mu = (w * y).sum() / w.sum()
        step = ((w ** 2 * (y - mu) ** 2).sum() - w.sum() + (w ** 2).sum() / w.sum()) / (w ** 2).sum()
        new = max(0.0, tau2 + step)
        if abs(new - tau2) < tol:
            tau2 = new
            break
        tau2 = new
    w = 1 / (v + tau2)
    mu = (w * y).sum() / w.sum()
    se = np.sqrt(1 / w.sum())
    return mu, se, tau2, mu / se


def load_qvals(qval_dir):
    """Gene-level q-values: a gs:// fdr_qvals dir (read with gcloud), or a local TSV(.gz) with gene_name, qval,
    cell_type (as written by cell_type_spec_runner.py's gather job)."""
    if not qval_dir.startswith('gs://'):
        return pd.read_csv(qval_dir, sep='\t')
    files = [f for f in sh(f"gcloud storage ls '{qval_dir}/'").decode().split() if f.endswith('_qval.tsv')]

    def read(f):
        d = pd.read_csv(io.BytesIO(sh(f"gcloud storage cat '{f}'")), sep='\t', usecols=['gene_name', 'qval'])
        return d.assign(cell_type=f.rsplit('/', 1)[1].removesuffix('_qval.tsv'))

    with ThreadPoolExecutor(16) as ex:
        return pd.concat(ex.map(read, files), ignore_index=True)


def load_expressed(path=None):
    """(gene, cell type) expressed in either cohort: the per-cell-type gene lists the eQTL scans were built from.
    path: a local TSV(.gz) with gene_name, cell_type (cell_type_spec_runner.py's gather job); default: read the lists
    from GCS with gcloud."""
    if path:
        d = pd.read_csv(path, sep='\t')
        return set(zip(d.gene_name, d.cell_type))

    def read(args):
        coh, ct = args
        raw = sh(f"gcloud storage cat '{GENE_LISTS.format(cohort=COHORT_DIRS[coh], ct=ct)}/*_gene_list.json'").decode()
        genes = set()
        for part in raw.replace('][', ']\n[').splitlines():
            if part.strip():
                genes.update(json.loads(part))
        return [(g, ct) for g in genes]

    with ThreadPoolExecutor(16) as ex:
        out = set()
        for chunk in ex.map(read, [(c, ct) for c in COHORTS for ct in CELL_TYPES]):
            out.update(chunk)
    return out


class SVLD:
    """Lazy per-chromosome SV dosages for both cohorts; signed r between federated SV pairs."""

    def __init__(self, cache_dir):
        import finemap_prep as fp  # noqa: PLC0415  (needs trtools)
        self.fp, self.cache_dir, self.loci = fp, cache_dir, {}
        os.makedirs(cache_dir, exist_ok=True)

    def _chrom(self, chrom):
        if chrom not in self.loci:
            self.loci.clear()
            per = {}
            for coh in COHORTS:
                path = os.path.join(self.cache_dir, f'{coh}_{chrom}.vcf.bgz')
                if not os.path.exists(path):
                    sh(f"gcloud storage cp -q '{SV_VCF.format(cohort=coh, chrom=chrom)}' '{path}'")
                loci = self.fp.load_sv_chrom(path)
                rows = np.arange(loci[0]['idx'].shape[0]) if loci else np.arange(0)
                per[coh] = {str(L['id']).upper(): self.fp.dosage(L, rows) for L in loci}
            self.loci[chrom] = per
        return self.loci[chrom]

    def r(self, chrom, sv_a, sv_b):
        """sv_a, sv_b: [bioheart_id, tob_id]. Returns (r_bioheart, r_tob, r_combined, r2_combined)."""
        per = self._chrom(chrom)
        rs, ns = [], []
        for k, coh in enumerate(COHORTS):
            a, b = per[coh].get(str(sv_a[k]).upper()), per[coh].get(str(sv_b[k]).upper())
            if a is None or b is None:
                rs.append(np.nan); ns.append(0); continue
            m = ~np.isnan(a) & ~np.isnan(b)
            if m.sum() < 20 or a[m].std() == 0 or b[m].std() == 0:
                rs.append(np.nan); ns.append(0); continue
            rs.append(float(np.corrcoef(a[m], b[m])[0, 1])); ns.append(int(m.sum()))
        w, r = np.array(ns, float), np.array(rs)
        if not w.sum():
            return rs[0], rs[1], np.nan, np.nan
        return rs[0], rs[1], float(np.nansum(r * w) / w.sum()), float(np.nansum(r ** 2 * w) / w.sum())


SS_COLS = ['gene_name', 'cell_type', 'chr', 'federate_id', 'pair_class', 'bioheart_vid', 'bioheart_pos', 'tob_vid',
           'coeff_meta_fixed', 'se_meta_fixed']


def main(sumstats, qval_dir, out_dir, rescue_effect='abs', vcf_cache=None, expressed_path=None):
    os.makedirs(out_dir, exist_ok=True)
    q = load_qvals(qval_dir)
    q = q[q.cell_type.isin(CELL_TYPES)]
    egene = set(zip(q.gene_name[q.qval < 0.05], q.cell_type[q.qval < 0.05]))
    expressed = load_expressed(expressed_path)
    print(f'{len(egene):,} eGene x cell type pairs; {len(expressed):,} expressed gene x cell type pairs')

    if sumstats.endswith('.parquet'):
        ss = pd.read_parquet(sumstats, columns=SS_COLS)
    else:   # TSV from the gather job: every test of every eGene's gene, all cell types
        ss = pd.read_csv(sumstats, sep='\t', usecols=SS_COLS, dtype={'chr': str})
    for k in ss.select_dtypes('category').columns:
        ss[k] = ss[k].astype(str)
    ss = ss[ss.cell_type.isin(CELL_TYPES)].copy()
    ss['z'] = ss.coeff_meta_fixed / ss.se_meta_fixed
    ss['absz'] = ss.z.abs()
    ss['sv'] = ss.bioheart_vid.str.upper() + '|' + ss.tob_vid.str.upper()
    tested = set(zip(ss.gene_name, ss.cell_type))
    print(f'{len(ss):,} SV x gene x cell type tests; {len(tested):,} gene x cell types tested; '
          f'tested but not on an expressed list: {len(tested - expressed):,}; '
          f'eGene x cell types not on an expressed list: {len(egene - expressed):,}')

    eg = ss[[k in egene for k in zip(ss.gene_name, ss.cell_type)]]
    mx = eg.groupby(['gene_name', 'cell_type']).absz.transform('max')
    leads = eg[eg.absz >= mx - TIE].copy()
    print(f'eQTLs (lead eSV rows): {len(leads):,} for {leads[["gene_name", "cell_type"]].drop_duplicates().shape[0]:,} '
          f'eGene x cell types; {leads.gene_name.nunique():,} eGenes')

    leads_by_gene = {g: d for g, d in leads.groupby('gene_name')}
    sub = ss[ss.gene_name.isin(leads_by_gene)].drop_duplicates(['gene_name', 'sv', 'cell_type'])
    sv_stats = dict(zip(zip(sub.gene_name, sub.sv, sub.cell_type), zip(sub.coeff_meta_fixed, sub.se_meta_fixed)))
    ld = SVLD(vcf_cache or os.path.join(out_dir, 'sv_vcfs'))
    ld_rows, ld_memo, rescue_rows = [], {}, []

    def corr(chrom, a, b):
        """(signed r, r2) between two SV pairs; r = 1 for the same SV."""
        if a == b:
            return 1.0, 1.0
        key = (a, b) if a <= b else (b, a)
        if key not in ld_memo:
            r1, r2_, rc, r2c = ld.r(chrom, key[0].split('|'), key[1].split('|'))
            ld_memo[key] = (rc, r2c)
            ld_rows.append(dict(chrom=chrom, sv_a=key[0], sv_b=key[1], r_bioheart=r1, r_tob=r2_, r=rc, r2=r2c))
        return ld_memo[key]

    results = []
    leads = leads.assign(chrom_=leads.chr.astype(str).str.split(',').str[0]).sort_values(['chrom_', 'gene_name', 'cell_type'])
    for i, row in enumerate(leads.itertuples(index=False)):
        if i % 5000 == 0:
            print(f'  {i:,} / {len(leads):,}', flush=True)
        g, A, sv_a, b_a, se_a, z_a = row.gene_name, row.cell_type, row.sv, row.coeff_meta_fixed, row.se_meta_fixed, row.z
        sign_a = np.sign(b_a)
        gl = leads_by_gene[g]
        scen, rescued = {}, []

        def put(ct, v):
            if v > scen.get(ct, -1):   # tied leads in B: keep the larger scenario (the eGene score is a max anyway)
                scen[ct] = v

        # eGene in B: 3 / 4 / 5 from the lead SV(s) in B
        for ct, d in gl[gl.cell_type != A].groupby('cell_type'):
            cs = [corr(row.chrom_, sv_a, s) for s in d.sv]
            r2s = np.array([c[1] if np.isfinite(c[1]) else -1 for c in cs])
            j = int(np.argmax(r2s))
            if r2s[j] < LD_R2:
                put(ct, 3)
            else:
                concordant = np.sign(d.coeff_meta_fixed.iloc[j]) * np.sign(cs[j][0]) == sign_a
                put(ct, 4 if concordant else 5)
        # not an eGene in B: 1 (not expressed) or 2 (expressed)
        eg_cts = set(gl.cell_type)
        for ct in CELL_TYPES:
            if ct != A and ct not in eg_cts:
                put(ct, 2 if (g, ct) in expressed else 1)
        # rescue 2 -> 4/5: A's lead SV in B, random-effects meta-analysis vs A's p
        lp_a = log10p(z_a)
        for ct in [k for k, v in scen.items() if v == 2]:
            if (g, sv_a, ct) not in sv_stats:
                continue
            b_b, se_b = sv_stats[(g, sv_a, ct)]
            y = np.abs([b_a, b_b]) if rescue_effect == 'abs' else np.array([b_a, b_b])
            mu, se_mu, tau2, z_m = rma_reml(y, [se_a ** 2, se_b ** 2])
            lp_m = log10p(z_m)
            hit = (lp_m - lp_a) < RESCUE_LOG10
            new = (4 if np.sign(b_b) == sign_a else 5) if hit else 2
            rescue_rows.append(dict(gene_name=g, cell_type=A, cell_type2=ct, sv=sv_a, coeff_main=b_a, se_main=se_a,
                                    log10p_main=lp_a, coeff_2=b_b, se_2=se_b, meta_estimate=mu, meta_se=se_mu, tau2=tau2,
                                    log10p_meta=lp_m, rescued=bool(hit), new_scenario=new))
            if hit:
                put(ct, new)
                rescued.append(ct)
        res = dict(gene_name=g, cell_type=A, chr=row.chrom_, federate_id=row.federate_id, pair_class=row.pair_class,
                   bioheart_vid=row.bioheart_vid, tob_vid=row.tob_vid, bioheart_pos=row.bioheart_pos,
                   coeff_meta_fixed=b_a, se_meta_fixed=se_a, z=z_a, log10p=lp_a, rescued_cts=','.join(rescued))
        res.update({f'scenario_{ct}': scen.get(ct, np.nan) for ct in CELL_TYPES})
        results.append(res)

    out = pd.DataFrame(results)
    sc = [f'scenario_{c}' for c in CELL_TYPES]
    out['specific_per_row'] = np.where((out[sc] == 4).any(axis=1), 'shared', 'specific')
    out['n_ct_shared'] = (out[sc] == 4).sum(axis=1) + 1                       # incl. the primary cell type
    gene_max = out.groupby('gene_name')[sc].max().max(axis=1).rename('max_scenario')
    out = out.merge(gene_max, left_on='gene_name', right_index=True)
    out['cell_spec_num'] = np.where(out.max_scenario.isin([1, 2, 3]), 1, out.n_ct_shared)
    out.to_csv(os.path.join(out_dir, 'esv_cell_type_specificity.tsv'), sep='\t', index=False)
    pd.DataFrame(ld_rows).to_csv(os.path.join(out_dir, 'sv_pair_ld.tsv'), sep='\t', index=False)
    pd.DataFrame(rescue_rows).to_csv(os.path.join(out_dir, 'rescue_tests.tsv'), sep='\t', index=False)
    rr = pd.DataFrame(rescue_rows)
    print(f'wrote {len(out):,} eQTL rows; {len(ld_rows):,} SV-SV LD pairs; {len(rr):,} scenario-2 meta-analyses, '
          f'{int(rr.rescued.sum()) if len(rr) else 0:,} rescued')


if __name__ == '__main__':
    import argparse  # noqa: PLC0415  (argparse, not click: the trtools image may not ship click)
    ap = argparse.ArgumentParser()
    ap.add_argument('--sumstats', required=True,
                    help='v2 full_sumstats.parquet, or the gather job TSV of every test of every eGene')
    ap.add_argument('--qvals', default=f'{V2}/fdr_qvals/using_acat',
                    help='gs:// fdr_qvals dir (read with gcloud) or a local TSV with gene_name, qval, cell_type')
    ap.add_argument('--expressed', default=None,
                    help='local TSV with gene_name, cell_type (default: read the gene lists from GCS with gcloud)')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--rescue-effect', default='abs', choices=['abs', 'signed'],
                    help='pool |beta| (default; rescues to 4 or 5) or signed beta (sensitivity)')
    ap.add_argument('--vcf-cache', default=None,
                    help='dir with {bioheart,tob}_{chrom}.vcf.bgz; missing files are copied from GCS with gcloud')
    a = ap.parse_args()
    main(a.sumstats, a.qvals, a.out_dir, a.rescue_effect, a.vcf_cache, a.expressed)
