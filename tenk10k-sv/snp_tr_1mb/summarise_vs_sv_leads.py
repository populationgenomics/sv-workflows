#!/usr/bin/env python3
"""
Summarise the SNP+TR 1-Mb re-test against the lead SV of every eGene x cell type pair.

For each meta file written by snp_tr_cis_runner.py ({meta-dir}/{chrom}_{chunk}.tsv.gz), one Hail Batch job keeps,
per gene x cell type x variant type (SNP / TR):
    - the top --top-n loci by meta-analysis |z| = |coeff_meta_fixed / se_meta_fixed|
      (|z| rather than p, because the strongest p-values underflow to 0);
    - the number of loci tested, and how many reach the lead SV's |z| (n_beats_sv).
The lead SV per gene x cell type is the max-|z| SV in the SV meta-analysis full_sumstats.parquet; the driver
computes it once and every job reads it.

Outputs, under the analysis output dir:
    sv_leads_all_pairs.parquet          lead SV per gene x cell type (every tested pair, not only significant)
    parts/{chrom}_{chunk}.top.parquet   per-file top loci (+ source_file, so a gene's chunk can be found later)
    parts/{chrom}_{chunk}.cnt.parquet   per-file counts
    snp_tr_top5.parquet, snp_tr_counts.parquet   all files combined (final job)
Existing parts are skipped, so a failed batch can be resubmitted.

analysis-runner --dataset tenk10k --access-level test --description "summarise SNP+TR 1Mb vs lead SVs" \
    --output-dir tenk10k-sv/sv/snp_tr_1mb/bioheart_n968_and_tob_n935/7pc/v1/summary_vs_sv_leads \
    python3 summarise_vs_sv_leads.py

For v2: v1 results + the v2 new-pair results, against the v2 lead SVs:

analysis-runner --dataset tenk10k --access-level test --description "summarise SNP+TR 1Mb vs v2 lead SVs" \
    --output-dir tenk10k-sv/sv/snp_tr_1mb/bioheart_n968_and_tob_n935/7pc/summary_vs_sv_leads_v2 \
    python3 summarise_vs_sv_leads.py \
    --meta-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/snp_tr_1mb/bioheart_n968_and_tob_n935/7pc/v1/meta,gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/snp_tr_1mb/bioheart_n968_and_tob_n935/7pc/v2_additions/meta \
    --sv-sumstats=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v2_F_filtered/full_sumstats.parquet
"""
import click

META_DIR = 'gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/snp_tr_1mb/bioheart_n968_and_tob_n935/7pc/v1/meta'
SV_SUMSTATS = ('gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/'
               'v1_F_filtered/full_sumstats.parquet')
COLS = ['gene_name', 'cell_type', 'variant_type', 'chr', 'pos', 'id', 'motif', 'coeff_1', 'se_1', 'coeff_2', 'se_2',
        'coeff_meta_fixed', 'se_meta_fixed', 'pval_meta_fixed', 'i2_meta', 'allele_frequency_1']
SV_COLS = ['gene_name', 'cell_type', 'federate_id', 'bioheart_vid', 'tob_vid', 'pair_class', 'chr', 'bioheart_pos',
           'coeff_1', 'se_1', 'coeff_2', 'se_2', 'coeff_meta_fixed', 'se_meta_fixed', 'pval_meta_fixed', 'i2_meta']
GROUP = ['gene_name', 'cell_type', 'variant_type']


def _blob(client, gs_path):
    bucket, _, name = gs_path.removeprefix('gs://').partition('/')
    return client.bucket(bucket).blob(name)


def _download(client, gs_path):
    import tempfile  # noqa: PLC0415
    fh = tempfile.NamedTemporaryFile(suffix='.' + gs_path.rsplit('.', 1)[-1], delete=False)
    _blob(client, gs_path).download_to_filename(fh.name)
    return fh.name


def _upload_parquet(client, df, gs_path):
    import tempfile  # noqa: PLC0415
    with tempfile.NamedTemporaryFile(suffix='.parquet') as fh:
        df.to_parquet(fh.name, index=False)
        _blob(client, gs_path).upload_from_filename(fh.name)


def sv_leads(sv_sumstats, out_path):
    """Max-|z| SV per gene x cell type from the SV meta-analysis."""
    import pandas as pd  # noqa: PLC0415
    from google.cloud import storage  # noqa: PLC0415
    client = storage.Client()
    sv = pd.read_parquet(_download(client, sv_sumstats), columns=SV_COLS)
    for k in sv.select_dtypes('category').columns:
        sv[k] = sv[k].astype(str)
    sv['absz'] = (sv.coeff_meta_fixed / sv.se_meta_fixed).abs()
    lead = sv.sort_values('absz', ascending=False).drop_duplicates(['gene_name', 'cell_type'])
    _upload_parquet(client, lead, out_path)
    print(f'{len(lead):,} gene x cell type pairs with a lead SV -> {out_path}')


def summarise_file(meta_path, leads_path, top_out, cnt_out, top_n, chunksize):
    """One meta file -> top loci and counts per gene x cell type x variant type."""
    import gzip  # noqa: PLC0415
    import os  # noqa: PLC0415
    import pandas as pd  # noqa: PLC0415
    from google.cloud import storage  # noqa: PLC0415

    client = storage.Client()
    if _blob(client, cnt_out).exists():
        print(f'{cnt_out} exists, skipping')
        return
    leads = pd.read_parquet(_download(client, leads_path), columns=['gene_name', 'cell_type', 'absz'])
    lead_z = leads.set_index(['gene_name', 'cell_type']).absz
    local = _download(client, meta_path)
    source = meta_path.removeprefix('gs://')   # full path, so the deep-dive can find a gene's file in any run
    tops, cnts = [], []
    with gzip.open(local, 'rt') as fh:
        for ch in pd.read_csv(fh, sep='\t', usecols=COLS, chunksize=chunksize, dtype={'motif': str, 'id': str}):
            ch['absz'] = (ch.coeff_meta_fixed / ch.se_meta_fixed).abs()
            sv_z = lead_z.reindex(pd.MultiIndex.from_frame(ch[['gene_name', 'cell_type']])).to_numpy()
            ch['beats_sv'] = ch.absz >= sv_z
            cnts.append(ch.groupby(GROUP).agg(n=('absz', 'size'), n_beats_sv=('beats_sv', 'sum')))
            tops.append(ch.drop(columns='beats_sv').sort_values('absz', ascending=False).groupby(GROUP).head(top_n))
    os.remove(local)
    # a gene x cell type never spans files, but it can span chunks of one file: re-aggregate
    cnt = pd.concat(cnts).groupby(level=[0, 1, 2]).sum().reset_index()
    top = pd.concat(tops).sort_values('absz', ascending=False).groupby(GROUP).head(top_n).assign(source_file=source)
    _upload_parquet(client, top, top_out)
    _upload_parquet(client, cnt, cnt_out)   # written last: its existence marks the file as done
    print(f'{source}: {len(cnt):,} gene x cell type x type groups')


def combine(part_paths, top_out, cnt_out):
    import pandas as pd  # noqa: PLC0415
    from google.cloud import storage  # noqa: PLC0415
    client = storage.Client()
    top = pd.concat(pd.read_parquet(_download(client, p + '.top.parquet')) for p in part_paths)
    cnt = pd.concat(pd.read_parquet(_download(client, p + '.cnt.parquet')) for p in part_paths)
    _upload_parquet(client, top, top_out)
    _upload_parquet(client, cnt, cnt_out)
    pairs = cnt[['gene_name', 'cell_type']].drop_duplicates()
    print(f'{len(part_paths)} files: {len(pairs):,} gene x cell type pairs, {len(top):,} top-locus rows')


@click.command()
@click.option('--meta-dir', default=META_DIR, help='one or more comma-separated meta dirs (e.g. v1 + v2_additions)')
@click.option('--sv-sumstats', default=SV_SUMSTATS)
@click.option('--top-n', default=5, type=int)
@click.option('--chunksize', default=1_000_000, type=int)
@click.option('--job-cpu', default=2, type=float)
def main(meta_dir, sv_sumstats, top_n, chunksize, job_cpu):
    from google.cloud import storage  # noqa: PLC0415
    from cpg_utils.hail_batch import get_batch, output_path  # noqa: PLC0415

    files = []
    for d in meta_dir.split(','):
        bucket, _, prefix = d.removeprefix('gs://').partition('/')
        found = sorted(f'gs://{bucket}/{b.name}' for b in storage.Client().list_blobs(bucket, prefix=prefix.rstrip('/') + '/')
                       if b.name.endswith('.tsv.gz'))
        print(f'{len(found)} meta files under {d}')
        files += found

    b = get_batch(name='summarise SNP+TR 1Mb vs lead SVs')
    leads_path = output_path('sv_leads_all_pairs.parquet', 'analysis')
    lj = b.new_python_job(name='sv_leads')
    lj.cpu(job_cpu); lj.memory('highmem')
    lj.call(sv_leads, sv_sumstats, leads_path)

    jobs, parts = [], []
    for f in files:
        # prefix with the run folder (e.g. v1 / v2_additions): chunk names repeat across runs
        name = f.rsplit('/', 3)[1] + '__' + f.rsplit('/', 1)[1].removesuffix('.tsv.gz')
        part = output_path(f'parts/{name}', 'analysis')
        parts.append(part)
        j = b.new_python_job(name=f'summarise_{name}')
        j.cpu(job_cpu); j.memory('highmem'); j.storage('10G')
        j.depends_on(lj)
        j.call(summarise_file, f, leads_path, part + '.top.parquet', part + '.cnt.parquet', top_n, chunksize)
        jobs.append(j)
    cj = b.new_python_job(name='combine')
    cj.memory('highmem'); cj.depends_on(*jobs)
    cj.call(combine, parts, output_path('snp_tr_top5.parquet', 'analysis'), output_path('snp_tr_counts.parquet', 'analysis'))
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
