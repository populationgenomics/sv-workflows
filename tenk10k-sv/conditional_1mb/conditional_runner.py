#!/usr/bin/env python3
# pylint: disable=too-many-locals,too-many-arguments,import-error,wrong-import-position
"""
Conditional analysis of SV-led eQTLs: for every eGene x cell type pair where the lead SV beats the best SNP and the
best TR in the gene's 1-Mb window, test every SNP, TR and SV in the window
    marginal          as is (the "before"),
    cond_lead_snptr   conditional on the lead SNP/TR,
    cond_lead_sv      conditional on the lead SV,
per cohort (BioHEART n = 968, TOB n = 935), then fixed-effect meta-analysis. See conditional_job.py.

Pair selection (same rule as the SV-led ranking, sv_led_1mb_ranking_v2.ipynb):
    eGene x cell type: gene-level q < 0.05 in --qval-dir (default: 7pc/v2_F_filtered)
    lead SV:           largest |z_meta| SV for the pair  ({summary}/sv_leads_all_pairs.parquet)
    best SNP / TR:     largest |z_meta| SNP / TR in the 1-Mb re-test ({summary}/snp_tr_top5.parquet)
    SV-led:            |z_SV| > |z_SNP| and |z_SV| > |z_TR|   (v2: 1,745 pairs, 537 genes)
    lead SNP/TR:       whichever of the best SNP and best TR has the larger |z|
--summary-dir is the output of ../snp_tr_1mb/summarise_vs_sv_leads.py for the same meta-analysis.

Outputs (per chromosome x chunk of --genes-per-job genes):
    {out}/assoc/{chrom}_{chunk:03d}.tsv.gz    every variant x condition: per-cohort and meta coeff/se/p, I2, and r2 of
                                              the variant with the lead SV and with the lead SNP/TR in each cohort
    {out}/summary/{chrom}_{chunk:03d}.tsv     one row per pair: z/p of the lead SV and lead SNP/TR under each
                                              condition, the top variant after each, LD between the two leads
    {out}/pairs.tsv                           the selected pairs and their conditioning variants
A chunk whose summary file exists is skipped, so a failed batch can be resubmitted. pairs.tsv is reused if present.
If the analysis-runner driver cannot read parquet, the first submission only runs a pair-selection job that writes
pairs.tsv; submit the same command again once it has finished.

Pilot (10 genes on chr22):
analysis-runner --dataset tenk10k --access-level test --description "SV-led conditional analysis pilot (v2)" \
    --output-dir tenk10k-sv/sv/conditional_1mb/bioheart_n968_and_tob_n935/7pc/v2_pilot \
    python3 conditional_runner.py --chromosomes=chr22 --max-genes=10

Full run:
analysis-runner --dataset tenk10k --access-level test --description "SV-led conditional analysis (v2)" \
    --output-dir tenk10k-sv/sv/conditional_1mb/bioheart_n968_and_tob_n935/7pc/v2 \
    python3 conditional_runner.py
"""
import json
import os
import sys

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.config import get_config

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'snp_tr_1mb'))
from snp_tr_cis_runner import SNP_VCF, SV_INPUTS, TR_VCF, cis_region  # noqa: E402

COHORTS = {'bioheart': 'bioheart_n968', 'tob': 'tob_n935'}
SV_VCF = ('gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/'
          '{cohort}_common_chr_specific/hail_filtered_{chrom}.vcf.bgz')
V2_META = 'gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v2_F_filtered'
V2_SUMMARY = ('gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/snp_tr_1mb/bioheart_n968_and_tob_n935/7pc/'
              'summary_vs_sv_leads_v2')
V2_META_SET = 'gs://cpg-tenk10k-test/sv/input_files/meta_analysis_set_v2_fixed_F_filtered.csv'


def sv_led_pairs(qval_dir, summary_dir, out_path=None):
    """SV-led eGene x cell type pairs with their lead SV and lead SNP/TR (the conditioning variants).

    Self-contained (imports inside, GCS via google-cloud-storage) so it can also run as a Hail Batch Python job when
    the driver cannot read parquet. Writes a TSV to out_path (gs:// or local) if given.
    """
    import io  # noqa: PLC0415
    from concurrent.futures import ThreadPoolExecutor  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415
    import pandas as pd  # noqa: PLC0415
    from google.cloud import storage  # noqa: PLC0415

    client = storage.Client()

    def blob(p):
        bucket, _, name = p.removeprefix('gs://').partition('/')
        return client.bucket(bucket).blob(name)

    def read_parquet(p, columns):
        d = pd.read_parquet(io.BytesIO(blob(p).download_as_bytes()), columns=columns)
        for k in d.select_dtypes('category').columns:
            d[k] = d[k].astype(str)
        return d

    bucket, _, prefix = qval_dir.removeprefix('gs://').partition('/')
    files = [b.name for b in client.list_blobs(bucket, prefix=prefix.rstrip('/') + '/') if b.name.endswith('_qval.tsv')]

    def read_q(name):
        d = pd.read_csv(io.BytesIO(client.bucket(bucket).blob(name).download_as_bytes()), sep='\t',
                        usecols=lambda c: c in ('gene_name', 'qval', 'chr'))
        return d.assign(cell_type=name.rsplit('/', 1)[1].removesuffix('_qval.tsv'))

    with ThreadPoolExecutor(16) as ex:
        q = pd.concat(ex.map(read_q, files), ignore_index=True)
    sig = q[q.qval < 0.05].copy()
    # genes with several lead SVs store a comma-separated chr list ('chr17,chr17'); all share one chromosome
    sig['chrom'] = sig.gene_name.map(sig.groupby('gene_name').chr.first().astype(str).str.split(',').str[0])
    lead = read_parquet(f'{summary_dir}/sv_leads_all_pairs.parquet',
                        ['gene_name', 'cell_type', 'federate_id', 'bioheart_vid', 'tob_vid', 'pair_class', 'absz'])
    top = read_parquet(f'{summary_dir}/snp_tr_top5.parquet', ['gene_name', 'cell_type', 'variant_type', 'id', 'pos', 'absz'])
    best = top.sort_values('absz', ascending=False).drop_duplicates(['gene_name', 'cell_type', 'variant_type'])
    w = best.pivot(index=['gene_name', 'cell_type'], columns='variant_type', values=['absz', 'id', 'pos'])
    w.columns = [f'{a}_{b}' for a, b in w.columns]
    p = (sig[['gene_name', 'cell_type', 'chrom']]
         .merge(lead.rename(columns={'absz': 'absz_SV'}), on=['gene_name', 'cell_type'], how='inner')
         .merge(w.reset_index(), on=['gene_name', 'cell_type'], how='inner'))
    for c in ('absz_SNP', 'absz_TR'):
        if c not in p:
            p[c] = np.nan
        p[c] = p[c].astype(float)
    # exact ties (identical genotypes, |z| equal to float precision) are not SV-led
    tol = 1e-4
    p['sv_led'] = (p.absz_SV > p.absz_SNP.fillna(-1) + tol) & (p.absz_SV > p.absz_TR.fillna(-1) + tol)
    p = p[p.sv_led].copy()
    use_snp = p.absz_SNP.fillna(-1) >= p.absz_TR.fillna(-1)
    p['lead_type'] = np.where(use_snp, 'SNP', 'TR')
    p['lead_id'] = np.where(use_snp, p.id_SNP, p.id_TR)
    p['lead_pos'] = np.where(use_snp, p.pos_SNP, p.pos_TR)
    p['lead_absz'] = np.where(use_snp, p.absz_SNP, p.absz_TR)
    p = p.sort_values(['chrom', 'gene_name', 'cell_type']).reset_index(drop=True)
    if out_path:
        if out_path.startswith('gs://'):
            blob(out_path).upload_from_string(p.to_csv(sep='\t', index=False), content_type='text/tab-separated-values')
        else:
            p.to_csv(out_path, sep='\t', index=False)
        print(f'{len(p):,} SV-led pairs -> {out_path}')
    return p


@click.command()
@click.option('--qval-dir', default=f'{V2_META}/fdr_qvals/using_acat')
@click.option('--summary-dir', default=V2_SUMMARY, help='summarise_vs_sv_leads.py output for the same meta-analysis')
@click.option('--meta-set', default=V2_META_SET, help='pairing table that meta-analysis pooled')
@click.option('--genes-per-job', default=10, type=int)
@click.option('--chromosomes', default=','.join(f'chr{i}' for i in range(1, 23)))
@click.option('--max-genes', default=0, type=int, help='pilot: cap genes per chromosome (0 = all)')
@click.option('--max-parallel-jobs', default=100, type=int)
@click.option('--storage', default='15G', help='chr1 TR VCF is 3.3 GB, SNP VCF 0.75 GB')
@click.option('--memory', default='standard')
@click.option('--write-pairs', default=None, help='only select the pairs, write them to this local TSV and exit')
def main(qval_dir, summary_dir, meta_set, genes_per_job, chromosomes, max_genes, max_parallel_jobs, storage, memory,
         write_pairs):
    if write_pairs:   # local, no Batch
        sv_led_pairs(qval_dir, summary_dir, write_pairs)
        return

    from cpg_utils.hail_batch import get_batch, output_path  # noqa: PLC0415
    pairs_path = output_path('pairs.tsv', 'analysis')
    if not to_path(pairs_path).exists():
        try:
            sv_led_pairs(qval_dir, summary_dir, pairs_path)
        except ImportError as e:   # driver image without a parquet engine: select in a Batch job, then rerun
            print(f'cannot read parquet here ({e}); submitting the pair-selection job. '
                  'Rerun this same command once it has finished.')
            b = get_batch(name='SV-led conditional analysis: select pairs')
            j = b.new_python_job(name='select SV-led pairs')
            j.memory('highmem')
            j.call(sv_led_pairs, qval_dir, summary_dir, pairs_path)
            b.run(wait=False)
            return
    with to_path(pairs_path).open() as fh:
        pairs = pd.read_csv(fh, sep='\t', dtype={'lead_id': str, 'id_SNP': str, 'id_TR': str})
    print(f'{len(pairs):,} SV-led eGene x cell type pairs, {pairs.gene_name.nunique():,} genes; '
          f'conditioning on a SNP in {(pairs.lead_type == "SNP").sum():,}, a TR in {(pairs.lead_type == "TR").sum():,}')

    srcs = {n: open(os.path.join(HERE, *p)).read() for n, p in {
        'fast_associatr.py': ('..', 'snp_tr_1mb', 'fast_associatr.py'),
        'meta_np.py': ('..', 'snp_tr_1mb', 'meta_np.py'),
        'finemap_prep.py': ('..', 'finemap_mcv', 'finemap_prep.py'),
        'conditional_job.py': ('conditional_job.py',)}.items()}

    b = get_batch(name='SV-led conditional analysis')
    meta_set_file = b.read_input(meta_set)
    jobs = []
    for chrom in chromosomes.split(','):
        genes = sorted(pairs.loc[pairs.chrom == chrom, 'gene_name'].unique())
        if max_genes:
            genes = genes[:max_genes]
        if not genes:
            continue
        snp = b.read_input_group(base=SNP_VCF.format(chrom=chrom), tbi=SNP_VCF.format(chrom=chrom) + '.tbi')
        tr = b.read_input_group(base=TR_VCF.format(chrom=chrom), tbi=TR_VCF.format(chrom=chrom) + '.tbi')
        svs = {c: b.read_input(SV_VCF.format(cohort=c, chrom=chrom)) for c in COHORTS}   # read end to end, no index
        for chunk, start in enumerate(range(0, len(genes), genes_per_job)):
            tag = f'{chrom}_{chunk:03d}'
            if to_path(output_path(f'summary/{tag}.tsv', 'analysis')).exists():
                continue
            manifest = {'chunk': tag, 'snp_vcf': str(snp.base), 'tr_vcf': str(tr.base), 'meta_set': str(meta_set_file),
                        'sv_vcf': {c: str(svs[c]) for c in COHORTS}, 'genes': []}
            for gene in genes[start:start + genes_per_job]:
                gp = pairs[pairs.gene_name == gene]
                entry = {'gene': gene, 'chrom': chrom, 'region': cis_region(gp.cell_type.iloc[0], chrom, gene), 'pairs': []}
                for r in gp.itertuples():
                    entry['pairs'].append({
                        'cell_type': r.cell_type,
                        'pheno': {c: str(b.read_input(
                            f'{SV_INPUTS}/{d}/7pc/pheno_cov_numpy/v1/{r.cell_type}/{chrom}/{gene}_pheno_cov.npy'))
                            for c, d in COHORTS.items()},
                        'lead_sv': {'bioheart_id': r.bioheart_vid, 'tob_id': r.tob_vid, 'pair_class': r.pair_class,
                                    'absz': float(r.absz_SV)},
                        'lead_nonsv': {'variant_type': r.lead_type, 'id': str(r.lead_id), 'pos': int(r.lead_pos),
                                       'absz': float(r.lead_absz)}})
                manifest['genes'].append(entry)

            j = b.new_job(name=f'conditional {tag}')
            j.image(get_config()['images']['trtools'])
            j.storage(storage)
            j.memory(memory)
            j.declare_resource_group(out={'assoc': '{root}.assoc.tsv.gz', 'summary': '{root}.summary.tsv'})
            j.command(
                'set -euo pipefail\nmkdir -p work && cd work\n'
                + ''.join(f"cat > {n} <<'__SRC_EOF__'\n{s}\n__SRC_EOF__\n" for n, s in srcs.items())
                + f"cat > manifest.json <<'__MAN_EOF__'\n{json.dumps(manifest)}\n__MAN_EOF__\n"
                + f"python3 conditional_job.py manifest.json {j.out['assoc']} {j.out['summary']} "
                  '2>&1 | { grep -v -i "warn" || true; }\n'
                + f"test -s {j.out['summary']}\n",
            )
            b.write_output(j.out['assoc'], output_path(f'assoc/{tag}.tsv.gz', 'analysis'))
            b.write_output(j.out['summary'], output_path(f'summary/{tag}.tsv', 'analysis'))   # marks the chunk done
            if len(jobs) >= max_parallel_jobs:
                j.depends_on(jobs[-max_parallel_jobs])
            jobs.append(j)
    print(f'{len(jobs)} jobs')
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
