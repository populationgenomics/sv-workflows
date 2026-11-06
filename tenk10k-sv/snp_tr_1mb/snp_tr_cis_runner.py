#!/usr/bin/env python3
# pylint: disable=too-many-locals,too-many-arguments
"""
Re-test SNPs and TRs in the SV analysis's 1-Mb cis windows, on the SV analysis's samples and covariates.

Scope: every eGene x cell type pair with gene-level q < 0.05 in the F-filtered SV meta-analysis (28,398 pairs,
7,612 genes). For each pair, every SNP (common_variant_snps) and TR (final-freeze tr_vcf) in the gene's 1-Mb cis
window is tested in BioHEART (n=968) and TOB (n=935) using exactly the SV run's inputs:

    cis window   gs://cpg-tenk10k-test/sv/input_files/bioheart_n968/7pc/cis_window_files/v1/{ct}/{chrom}/{gene}_1000000bp.bed
    pheno_cov    gs://cpg-tenk10k-test/sv/input_files/{cohort}/7pc/pheno_cov_numpy/v1/{ct}/{chrom}/{gene}_pheno_cov.npy

so SNP/TR and SV betas, SEs and p-values come from the same individuals, expression values, covariates and window.

Engine: fast_associatr.py, an associaTR 6.0.1-equivalent (see README.md for the equivalence test), followed by
meta_np.metagen_dl, a numpy re-implementation of meta::metagen(..., random = TRUE, method.tau = "DL") as used
by meta_runner_sv.py. Both VCFs contain both cohorts, so one genotype pass per gene serves every cell type in
both cohorts.

Jobs: one per (chromosome, chunk of --genes-per-job genes). Each job localises the chromosome's SNP and TR VCFs
once, then loops over its genes. Outputs per job (gzipped TSV):

    {out}/meta/{chrom}_{chunk:03d}.tsv.gz       one row per gene x cell type x locus tested in both cohorts:
                                                  gene_name, cell_type, variant_type (SNP/TR), chr, pos, id, motif,
                                                  period, ref_len, n_samples_tested_1/2, coeff_1/2, se_1/2, pval_1/2,
                                                  r2_1/2, allele_frequency_1/2, + the meta_runner_sv.py meta columns
    {out}/per_cohort/{chrom}_{chunk:03d}.tsv.gz associaTR-format rows for every locus (incl. filtered), per cohort
Cohort suffix: _1 = BioHEART, _2 = TOB (as in the SV meta-analysis).

analysis-runner --dataset tenk10k --access-level test --description "SNP+TR 1Mb cis re-test on SV samples" \
    --output-dir tenk10k-sv/sv/snp_tr_1mb/bioheart_n968_and_tob_n935/7pc/v1 \
    python3 snp_tr_cis_runner.py --genes-per-job=50

Add --chromosomes=chr22 --max-genes=20 for a pilot.
"""
import json
import os
from concurrent.futures import ThreadPoolExecutor

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.config import get_config

HERE = os.path.dirname(os.path.abspath(__file__))
SV_BASE = 'gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v1_F_filtered'
SV_INPUTS = 'gs://cpg-tenk10k-test/sv/input_files'
COHORTS = {'bioheart': ('bioheart_n968', '1'), 'tob': ('tob_n935', '2')}
SNP_VCF = 'gs://cpg-bioheart-test/tenk10k/str/associatr/common_variant_snps/hail_filtered_{chrom}.vcf.bgz'
TR_VCF = 'gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific/hail_filtered_{chrom}.vcf.bgz'

# Runs inside each job (trtools image). Reads manifest.json, writes the two gzipped TSVs.
JOB_SCRIPT = r'''
import gzip, json, os, sys
import numpy as np
import pandas as pd
import fast_associatr as fa
from meta_np import metagen_dl

# input paths are written as ${BATCH_TMPDIR}/... inside a quoted heredoc, so expand them here
man = json.loads(os.path.expandvars(open(sys.argv[1]).read()))
meta_out, coh_out = gzip.open(sys.argv[2], 'wt'), gzip.open(sys.argv[3], 'wt')
meta_header = coh_header = True
COLS = ['chrom', 'pos', 'alleles', 'n_samples_tested', 'locus_filtered', 'p', 'coeff', 'se', 'r2',
        'motif', 'period', 'ref_len', 'allele_frequency']
for g in man['genes']:
    for vtype, vcf in (('SNP', man['snp_vcf']), ('TR', man['tr_vcf'])):
        phenos = {f"{coh}|{ct}": path for ct, by_coh in g['pheno'].items() for coh, path in by_coh.items()}
        res = fa.run_window(vcf, g['region'], 'eh', phenos)
        ids = fa.LAST_IDS
        for ct in g['pheno']:
            per = {}
            for coh, sfx in (('bioheart', '1'), ('tob', '2')):
                d = pd.DataFrame(res[f'{coh}|{ct}'], columns=COLS)
                d.insert(0, 'id', ids)
                per[sfx] = d
                o = d.assign(gene_name=g['gene'], cell_type=ct, cohort=coh, variant_type=vtype)
                o.to_csv(coh_out, sep='\t', index=False, header=coh_header); coh_header = False
            a, b = per['1'], per['2']
            ok = (a.locus_filtered.astype(str) == 'False') & (b.locus_filtered.astype(str) == 'False')
            if not ok.any():
                continue
            m = metagen_dl(np.column_stack([a.coeff[ok], b.coeff[ok]]), np.column_stack([a.se[ok], b.se[ok]]))
            out = pd.DataFrame({
                'gene_name': g['gene'], 'cell_type': ct, 'variant_type': vtype, 'chr': a.chrom[ok], 'pos': a.pos[ok],
                'id': a.id[ok], 'motif': a.motif[ok], 'period': a.period[ok], 'ref_len': a.ref_len[ok],
                'n_samples_tested_1': a.n_samples_tested[ok], 'n_samples_tested_2': b.n_samples_tested[ok],
                'coeff_1': a.coeff[ok], 'se_1': a.se[ok], 'pval_1': a.p[ok].astype(float), 'r2_1': a.r2[ok],
                'coeff_2': b.coeff[ok], 'se_2': b.se[ok], 'pval_2': b.p[ok].astype(float), 'r2_2': b.r2[ok],
                'allele_frequency_1': a.allele_frequency[ok], 'allele_frequency_2': b.allele_frequency[ok], **m})
            out.to_csv(meta_out, sep='\t', index=False, header=meta_header); meta_header = False
    print(f"{g['gene']} done ({len(g['pheno'])} cell types)", flush=True)
meta_out.close(); coh_out.close()
'''


def egene_pairs(qval_dir):
    """eGene x cell type pairs at gene-level q < 0.05, with chromosome."""
    files = [f for f in to_path(qval_dir).glob('*_qval.tsv')]

    def read(f):
        with f.open() as h:
            d = pd.read_csv(h, sep='\t', usecols=lambda c: c in ('gene_name', 'qval', 'chr', 'chrom'))
        return d.assign(cell_type=f.name.removesuffix('_qval.tsv'))

    with ThreadPoolExecutor(16) as ex:
        q = pd.concat(ex.map(read, files), ignore_index=True)
    return q[q.qval < 0.05]


def gene_chrom(pairs):
    """Chromosome per gene from the SV pheno_cov layout (the q-value files may not carry it)."""
    if 'chr' in pairs or 'chrom' in pairs:
        col = 'chr' if 'chr' in pairs else 'chrom'
        return pairs.groupby('gene_name')[col].first().to_dict()
    out = {}
    for ct in sorted(pairs.cell_type.unique()):
        root = f'{SV_INPUTS}/bioheart_n968/7pc/pheno_cov_numpy/v1/{ct}'
        need = set(pairs.loc[pairs.cell_type == ct, 'gene_name']) - set(out)
        if not need:
            continue
        for chrom_dir in to_path(root).iterdir():
            for f in chrom_dir.glob('*_pheno_cov.npy'):
                gname = f.name.removesuffix('_pheno_cov.npy')
                if gname in need:
                    out[gname] = chrom_dir.name
    return out


def cis_region(ct, chrom, gene):
    p = f'{SV_INPUTS}/bioheart_n968/7pc/cis_window_files/v1/{ct}/{chrom}/{gene}_1000000bp.bed'
    with to_path(p).open() as h:
        c, s, e = h.readline().split('\t')[:3]
    return f'{c}:{int(s)}-{int(e)}'


@click.command()
@click.option('--qval-dir', default=f'{SV_BASE}/fdr_qvals/using_acat')
@click.option('--genes-per-job', default=50, type=int)
@click.option('--chromosomes', default=','.join(f'chr{i}' for i in range(1, 23)))
@click.option('--max-genes', default=0, type=int, help='pilot: cap genes per chromosome (0 = all)')
@click.option('--max-parallel-jobs', default=200, type=int)
@click.option('--storage', default='12G', help='job disk: chr1 TR VCF is 3.3 GB, SNP VCF 0.75 GB')
@click.option('--memory', default='standard')
def main(qval_dir, genes_per_job, chromosomes, max_genes, max_parallel_jobs, storage, memory):
    from cpg_utils.hail_batch import get_batch, output_path

    pairs = egene_pairs(qval_dir)
    print(f'{len(pairs):,} eGene x cell type pairs, {pairs.gene_name.nunique():,} genes')
    chrom_of = gene_chrom(pairs)
    pairs['chrom'] = pairs.gene_name.map(chrom_of)
    assert pairs.chrom.notna().all(), pairs[pairs.chrom.isna()].head()

    fast_src = open(os.path.join(HERE, 'fast_associatr.py')).read()
    meta_src = open(os.path.join(HERE, 'meta_np.py')).read()

    b = get_batch(name='SNP+TR 1Mb cis re-test on SV samples')
    jobs = []
    for chrom in chromosomes.split(','):
        genes = sorted(pairs.loc[pairs.chrom == chrom, 'gene_name'].unique())
        if max_genes:
            genes = genes[:max_genes]
        if not genes:
            continue
        snp = b.read_input_group(base=SNP_VCF.format(chrom=chrom), tbi=SNP_VCF.format(chrom=chrom) + '.tbi')
        tr = b.read_input_group(base=TR_VCF.format(chrom=chrom), tbi=TR_VCF.format(chrom=chrom) + '.tbi')
        for chunk, start in enumerate(range(0, len(genes), genes_per_job)):
            meta_path = output_path(f'meta/{chrom}_{chunk:03d}.tsv.gz', 'analysis')
            if to_path(meta_path).exists():
                continue
            batch_genes = genes[start:start + genes_per_job]
            j = b.new_job(name=f'snp_tr_1mb_{chrom}_{chunk:03d}')
            j.image(get_config()['images']['trtools'])
            j.storage(storage)
            j.memory(memory)
            manifest = {'snp_vcf': str(snp.base), 'tr_vcf': str(tr.base), 'genes': []}
            for gene in batch_genes:
                cts = sorted(pairs.loc[(pairs.gene_name == gene), 'cell_type'])
                entry = {'gene': gene, 'region': cis_region(cts[0], chrom, gene), 'pheno': {}}
                for ct in cts:
                    entry['pheno'][ct] = {
                        coh: str(b.read_input(f'{SV_INPUTS}/{d}/7pc/pheno_cov_numpy/v1/{ct}/{chrom}/{gene}_pheno_cov.npy'))
                        for coh, (d, _) in COHORTS.items()}
                manifest['genes'].append(entry)
            j.declare_resource_group(out={'meta': '{root}.meta.tsv.gz', 'cohort': '{root}.cohort.tsv.gz'})
            j.command(
                'set -euo pipefail\nmkdir -p work && cd work\n'
                f"cat > fast_associatr.py <<'__FAST_EOF__'\n{fast_src}\n__FAST_EOF__\n"
                f"cat > meta_np.py <<'__META_EOF__'\n{meta_src}\n__META_EOF__\n"
                f"cat > job.py <<'__JOB_EOF__'\n{JOB_SCRIPT}\n__JOB_EOF__\n"
                f"cat > manifest.json <<'__MAN_EOF__'\n{json.dumps(manifest)}\n__MAN_EOF__\n"
                f"python3 job.py manifest.json {j.out['meta']} {j.out['cohort']}\n",
            )
            b.write_output(j.out['meta'], meta_path)
            b.write_output(j.out['cohort'], output_path(f'per_cohort/{chrom}_{chunk:03d}.tsv.gz', 'analysis'))
            if len(jobs) >= max_parallel_jobs:
                j.depends_on(jobs[-max_parallel_jobs])
            jobs.append(j)
    print(f'{len(jobs)} jobs')
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
