"""Writes conditional_plots.ipynb: summaries and before/after association plots for conditional_runner.py output.

python3 make_conditional_plots_nb.py            # then open conditional_plots.ipynb (VM or local)
"""
import nbformat as nbf

md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
C = [md(r'''# SV-led eQTLs: conditional analysis (v2)

For every eGene × cell type pair where the lead SV has a larger |z| than the best SNP and the best TR in the gene's
1-Mb window (1,735 pairs), every SNP, TR and SV in the window was re-tested three times (`conditional_runner.py`):

| condition | model, per cohort, then fixed-effect meta-analysis |
|---|---|
| `marginal` | expression ~ variant + covariates (the "before"; reproduces the eQTL scans) |
| `cond_lead_snptr` | ... + lead SNP/TR dosage |
| `cond_lead_sv` | ... + lead SV dosage |

**Reading the two key tests**
- **SV | lead SNP/TR:** does the SV keep its association once the strongest SNP/TR is in the model?
- **lead SNP/TR | SV:** does the SNP/TR signal disappear once the SV is in the model?

| SV \| SNP/TR | SNP/TR \| SV | call |
|---|---|---|
| significant | not significant | **SV explains** the SNP/TR signal |
| significant | significant | **independent** signals (both needed) |
| not significant | not significant | **not separable** (proxies of each other, high LD) |
| not significant | significant | SNP/TR explains the SV signal |

"Significant" is p < `P_SIG` (default 1e-5); the table is also shown at 5e-8 and 0.05. p-values are computed from
z = β/SE (two-sided normal), because the strongest p-values underflow to 0.'''),
code(r'''import glob, gzip, io, os, re, subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats

RESULTS = 'gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/conditional_1mb/bioheart_n968_and_tob_n935/7pc/v2'
P_SIG = 1e-5
GTF = os.path.expanduser('~/sv/exonic_deletions/gencode.v44.basic.annotation.gtf.gz')   # optional gene track
GENE_ANN = 'gene_ann.csv'                                                                  # optional symbols
OUT = 'conditional_plots'; os.makedirs(OUT, exist_ok=True)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_colwidth', 60)
mpl.rcParams.update({'font.family': 'sans-serif', 'font.size': 6, 'axes.titlesize': 6, 'axes.labelsize': 6,
                     'xtick.labelsize': 5, 'ytick.labelsize': 5, 'legend.fontsize': 5, 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False})
CONDS = ['marginal', 'cond_lead_snptr', 'cond_lead_sv']


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, executable='/bin/bash').stdout


def ls(pattern):
    if RESULTS.startswith('gs://'):
        return sorted(sh(f"gcloud storage ls '{RESULTS}/{pattern}'").decode().split())
    return sorted(glob.glob(f'{RESULTS}/{pattern}'))


def read_bytes(path):
    if path.startswith('gs://'):
        return sh(f"gcloud storage cat '{path}'")
    return open(path, 'rb').read()


def mlogp(z):
    """-log10 two-sided p from z, without underflow."""
    return -(stats.norm.logsf(np.abs(z)) + np.log(2)) / np.log(10)'''),
md('## 1. Load the per-pair summaries'),
code(r'''S = pd.concat([pd.read_csv(io.BytesIO(read_bytes(f)), sep='\t') for f in ls('summary/*.tsv')], ignore_index=True)
skipped = S[S.note.notna() & (S.note.astype(str) != '')]
S = S[~S.index.isin(skipped.index)].copy()
print(f'{len(S):,} pairs analysed ({S.gene_name.nunique():,} genes); skipped: {len(skipped)}')
if len(skipped):
    print(skipped[['gene_name', 'cell_type', 'note']].to_string(index=False))
S['svtype'] = S.lead_sv_bioheart.str.extract(r'final_cleanup_([A-Z]+)_', expand=False)
S['r2_sv_lead'] = S[['r2_sv_lead_1', 'r2_sv_lead_2']].mean(axis=1)
for c in CONDS:
    S[f'mlp_sv_{c}'] = mlogp(S[f'z_sv_{c}'])
    S[f'mlp_lead_{c}'] = mlogp(S[f'z_lead_{c}'])
if os.path.exists(GENE_ANN):
    S = S.merge(pd.read_csv(GENE_ANN)[['gene_name', 'symbol']], on='gene_name', how='left')
else:
    S['symbol'] = np.nan
S['symbol'] = S.symbol.fillna(S.gene_name)
# sanity: the marginal pass reproduces the selection inputs
d = (S.z_sv_marginal.abs() - S.input_absz_sv).abs().max(), (S.z_lead_marginal.abs() - S.input_absz_nonsv).abs().max()
print(f'max |marginal |z| - input |z||: SV {d[0]:.2e}, lead SNP/TR {d[1]:.2e}')'''),
md('## 2. Classify each pair'),
code(r'''def classify(S, p_sig):
    t = -np.log10(p_sig)
    sv_stays = S.mlp_sv_cond_lead_snptr > t
    lead_stays = S.mlp_lead_cond_lead_sv > t
    return np.select([sv_stays & ~lead_stays, sv_stays & lead_stays, ~sv_stays & ~lead_stays],
                     ['SV explains', 'independent', 'not separable'], 'SNP/TR explains')


ORDER = ['SV explains', 'independent', 'not separable', 'SNP/TR explains']
S['call'] = classify(S, P_SIG)
tab = pd.DataFrame({f'p < {p:g}': pd.Series(classify(S, p)).value_counts().reindex(ORDER, fill_value=0)
                    for p in (5e-8, P_SIG, 0.05)})
print(f'{len(S):,} SV-led pairs; rows = call, columns = significance threshold used after conditioning')
print(tab.to_string())
print('\nby lead SNP/TR type (p < P_SIG):')
print(pd.crosstab(S.lead_nonsv_type, S.call).reindex(columns=ORDER, fill_value=0).to_string())
print('\nby SV type (p < P_SIG):')
print(pd.crosstab(S.svtype, S.call).reindex(columns=ORDER, fill_value=0).to_string())
S['r2_bin'] = pd.cut(S.r2_sv_lead, [-0.01, 0.2, 0.5, 0.8, 1.0], labels=['<0.2', '0.2-0.5', '0.5-0.8', '>=0.8'])
print('\nby LD between the lead SV and the lead SNP/TR (r2, mean of cohorts):')
print(pd.crosstab(S.r2_bin, S.call).reindex(columns=ORDER, fill_value=0).to_string())
genes = S.groupby('gene_name').call.agg(lambda c: 'SV explains' if (c == 'SV explains').any() else c.iloc[0])
print(f'\neGenes with >= 1 cell type where the SV explains the SNP/TR signal: {(genes == "SV explains").sum():,} of {len(genes):,}')
S.to_csv(f'{OUT}/conditional_summary_classified.tsv', sep='\t', index=False)'''),
md(r'''## 3. Before vs after, all pairs
Left: the SV, before and after conditioning on the lead SNP/TR. Right: the lead SNP/TR, before and after
conditioning on the SV. Points below the dashed line lost significance (p < `P_SIG`).'''),
code(r'''COL = {'SV explains': '#d7191c', 'independent': '#2c7bb6', 'not separable': '#7f7f7f', 'SNP/TR explains': '#fdae61'}
fig, axes = plt.subplots(1, 2, figsize=(150 / 25.4, 68 / 25.4))
t = -np.log10(P_SIG)
for ax, (x, y, lab) in zip(axes, [('mlp_sv_marginal', 'mlp_sv_cond_lead_snptr', 'lead SV'),
                                  ('mlp_lead_marginal', 'mlp_lead_cond_lead_sv', 'lead SNP/TR')]):
    for c in ORDER:
        d = S[S.call == c]
        ax.scatter(d[x], d[y], s=4, lw=0, alpha=0.7, color=COL[c], label=f'{c} ({len(d):,})')
    lim = max(S[x].max(), S[y].max()) * 1.05
    ax.plot([0, lim], [0, lim], color='0.6', lw=0.5)
    ax.axhline(t, color='0.3', lw=0.5, ls='--')
    ax.set_xscale('symlog', linthresh=10); ax.set_yscale('symlog', linthresh=10)
    ax.set_xlabel(f'{lab}: −log10 p, marginal')
    other = 'lead SNP/TR' if lab == 'lead SV' else 'lead SV'
    ax.set_ylabel(f'{lab}: −log10 p, conditional on {other}')
axes[0].legend(frameon=False, loc='upper left')
fig.tight_layout(); fig.savefig(f'{OUT}/before_after_all_pairs.pdf', bbox_inches='tight'); plt.show()'''),
md(r'''## 4. Association plots before and after conditioning

Three panels per pair, same x and y axes: marginal; conditional on the lead SNP/TR; conditional on the lead SV.
SNPs are coloured by r² with the lead SV (mean of the two cohorts); TRs are grey triangles; other SVs are purple
squares; the lead SV is the red diamond and the lead SNP/TR is outlined in black. The conditioning variant is absent
from its own conditional panel (its position is marked by a dotted line).'''),
code(r'''LD_BINS = [(0.8, '#d7191c'), (0.6, '#fdae61'), (0.4, '#a6d96a'), (0.2, '#1a9641'), (0.0, '#2b83ba')]
ld_col = lambda r2: '#bdbdbd' if not np.isfinite(r2) else next(col for cut, col in LD_BINS if r2 >= cut)
_cache = {}


def pair_rows(gene, ct):
    chunk = S.loc[(S.gene_name == gene) & (S.cell_type == ct), 'chunk'].iloc[0]
    if chunk not in _cache:
        _cache.clear()
        _cache[chunk] = pd.read_csv(io.BytesIO(read_bytes(f'{RESULTS}/assoc/{chunk}.tsv.gz')), sep='\t',
                                    compression='gzip', dtype={'id': str, 'tob_id': str})
    d = _cache[chunk]
    d = d[(d.gene_name == gene) & (d.cell_type == ct)].copy()
    d['mlp'] = mlogp(d.coeff_meta_fixed / d.se_meta_fixed)
    d['r2'] = d[['r2_sv_1', 'r2_sv_2']].mean(axis=1)
    return d


def gene_track(gene_id):
    if not os.path.exists(GTF):
        return None
    out = sh(f"zcat {GTF} | awk -F'\\t' '$3==\"gene\" || $3==\"transcript\"' | grep '{gene_id}' || true").decode()
    g = mane = None
    for line in out.splitlines():
        f = line.split('\t')
        if f[2] == 'gene':
            g = (int(f[3]), int(f[4]), 1 if f[6] == '+' else -1)
        elif 'MANE_Select' in f[8]:
            mane = (int(f[3]), int(f[4]))
    return (g, mane or (g[:2] if g else None)) if g else None


def plot_pair(gene, ct, save=True):
    s = S[(S.gene_name == gene) & (S.cell_type == ct)].iloc[0]
    d = pair_rows(gene, ct)
    lead_sv, lt, lid = s.lead_sv_bioheart, s.lead_nonsv_type, str(s.lead_nonsv_id)
    lead_pos = d[(d.variant_type == lt) & (d.id == lid)].pos
    lead_pos = lead_pos.iloc[0] if len(lead_pos) else np.nan
    ymax = d.mlp.max() * 1.12
    track = gene_track(gene)
    fig, axes = plt.subplots(3, 1, figsize=(150 / 25.4, 120 / 25.4), sharex=True)
    titles = {'marginal': 'marginal',
              'cond_lead_snptr': f'conditional on lead {lt} {lid}',
              'cond_lead_sv': 'conditional on lead SV'}
    for ax, cond in zip(axes, CONDS):
        x = d[d.condition == cond]
        snp, tr, sv = x[x.variant_type == 'SNP'], x[x.variant_type == 'TR'], x[x.variant_type == 'SV']
        ax.scatter(tr.pos / 1e6, tr.mlp, s=6, marker='^', facecolor='none',
                   edgecolor=np.where(tr.mlp < 2, '0.82', '0.45'), lw=0.35, zorder=1)
        snp = snp.assign(_o=snp.r2.fillna(-1)).sort_values(['_o', 'mlp'])
        ax.scatter(snp.pos / 1e6, snp.mlp, s=np.where(snp.mlp < 2, 1.5, 4), c=[ld_col(r) for r in snp.r2], lw=0, zorder=2)
        oth = sv[sv.id != lead_sv]
        ax.scatter(oth.pos / 1e6, oth.mlp, s=9, marker='s', facecolor='none', edgecolor='#6a3d9a', lw=0.5, zorder=3)
        ls_ = sv[sv.id == lead_sv]
        if len(ls_):
            ax.scatter(ls_.pos / 1e6, ls_.mlp, s=34, marker='D', color='#e31a1c', edgecolor='k', lw=0.5, zorder=5)
        ln = x[(x.variant_type == lt) & (x.id == lid)]
        if len(ln):
            ax.scatter(ln.pos / 1e6, ln.mlp, s=30, marker='o' if lt == 'SNP' else '^', color=ld_col(ln.r2.iloc[0]),
                       edgecolor='k', lw=0.9, zorder=6)
        if cond == 'cond_lead_snptr' and np.isfinite(lead_pos):
            ax.axvline(lead_pos / 1e6, color='0.4', lw=0.5, ls=':')
        if cond == 'cond_lead_sv':
            ax.axvline(d[(d.id == lead_sv)].pos.iloc[0] / 1e6, color='#e31a1c', lw=0.5, ls=':')
        ax.axhline(-np.log10(P_SIG), color='0.6', lw=0.4, ls='--')
        ax.set_ylim(-0.12 * ymax if track else 0, ymax)
        ax.set_ylabel('−log10 p')
        z_sv = s[f'z_sv_{cond}']; z_ld = s[f'z_lead_{cond}']
        bits = []
        if np.isfinite(z_sv): bits.append(f'SV p = {10 ** -mlogp(z_sv):.1e}')
        if np.isfinite(z_ld): bits.append(f'lead {lt} p = {10 ** -mlogp(z_ld):.1e}')
        ax.set_title(f"{titles[cond]}   ({'; '.join(bits)})", loc='left')
        if track:
            (g0, g1, strand), (m0, m1) = track
            yy = -0.06 * ymax
            ax.plot([g0 / 1e6, g1 / 1e6], [yy] * 2, color='k', lw=0.6, solid_capstyle='butt')
            ax.plot([m0 / 1e6, m1 / 1e6], [yy] * 2, color='k', lw=2.2, solid_capstyle='butt')
            ax.set_yticks([v for v in ax.get_yticks() if 0 <= v <= ymax])
    axes[-1].set_xlabel(f'{s.chrom} position (Mb)')
    fig.suptitle(f"{s.symbol} · {ct.replace('_', ' ')} · {s.call} (r² SV–lead {lt} = {s.r2_sv_lead:.2f})",
                 x=0.01, ha='left', fontsize=7, y=0.995)
    h = [Line2D([], [], marker='o', ls='', color=col, ms=3, label=f'SNP r² with SV ≥ {cut}') for cut, col in LD_BINS]
    h += [Line2D([], [], marker='^', ls='', mfc='none', mec='0.35', ms=3, label='TR'),
          Line2D([], [], marker='s', ls='', mfc='none', mec='#6a3d9a', ms=3, label='other SVs'),
          Line2D([], [], marker='D', ls='', color='#e31a1c', mec='k', ms=4, label='lead SV'),
          Line2D([], [], marker='o', ls='', color='white', mec='k', ms=4, label=f'lead SNP/TR')]
    axes[0].legend(handles=h, frameon=False, loc='upper left', bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    if save:
        fig.savefig(f'{OUT}/{s.symbol}_{ct}_conditional.pdf', bbox_inches='tight')
    plt.show()'''),
md(r'''### Examples: the strongest pairs in each class
Edit `PICK` to plot any pair, e.g. `PICK = [('ENSG00000155530', 'NK')]`.'''),
code(r'''PICK = []
for c in ORDER[:3]:
    top = S[S.call == c].sort_values('mlp_sv_marginal', ascending=False).drop_duplicates('gene_name').head(3)
    PICK += list(zip(top.gene_name, top.cell_type))
print(S.set_index(['gene_name', 'cell_type']).loc[PICK][['symbol', 'svtype', 'call', 'lead_nonsv_type', 'r2_sv_lead',
      'z_sv_marginal', 'z_sv_cond_lead_snptr', 'z_lead_marginal', 'z_lead_cond_lead_sv']].round(2).to_string())
for gene, ct in PICK:
    plot_pair(gene, ct)'''),
]
nb = nbf.v4.new_notebook(); nb.cells = C
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
nbf.write(nb, 'conditional_plots.ipynb')
print('written conditional_plots.ipynb')
