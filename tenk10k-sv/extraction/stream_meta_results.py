#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
Stream the full meta-analysis results straight out of GCS and distil them to one compact Parquet.

Reads
    gs://.../meta_analysis/{cohorts}/{npc}/{version}/meta_results/{cell_type}/{chrom}/{gene}_*.tsv

and never writes the source TSVs to disk. Under analysis-runner it fans out one Hail Batch job per
cell type; each job reads that cell type's per-gene files in parallel with google-cloud-storage,
projects them to the columns that matter and writes one part file

    {parts-dir}/{cell_type}.parquet

and a final job, which waits for all of them, concatenates the parts into --out (one or more row
groups per cell type, so `cell_type` filters still skip most of the file). Part files that already
exist are not recomputed, so a failed batch can be resubmitted and only redoes what is missing.
--local runs the same two steps in-process, one cell type at a time.

WHY THIS EXISTS
---------------
`extract_lead_sv_sumstats.py` exports only the LEAD SV per (cell type, gene). That is the right
size for sharing, but it cannot answer anything that needs the tested universe:

  - enrichment tests need the tested set as the denominator and to build a matched null;
  - every effect-size number computed from lead SVs is conditioned on being the lead, which biases
    |beta| upward;
  - QQ plots and lambda_GC need every tested p-value, not the minimum per gene.

The tested universe is NOT derivable from SV intervals and cis-window geometry: interval overlap
predicts ~19.8 testable SVs per gene against ~12.1 actually tested (median ratio 0.63, 5th-95th
percentile 0.31-0.90), because federation between the two cohorts and associaTR's per-cell-type
filters both drop loci. Hence reading the real thing.

SIZE
----
~379k files, ~9.2 rows each, ~3.5M rows total. One cell type takes ~1-3 min at 64 threads, so the
batch finishes in a few minutes plus the combine step, and the output is tens of MB -- small enough
that every later query is local and instant.

PRIVACY
-------
Every retained column is a regression summary statistic, a variant identifier, an allele
frequency or a sample COUNT. No per-sample value is read or written.

analysis-runner --dataset "tenk10k" --description "stream SV meta results to parquet" --access-level "test" \
    --output-dir "tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v1_F_filtered" \
    python3 stream_meta_results.py \
    --meta-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v1_F_filtered/meta_results \
    --out=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v1_F_filtered/full_sumstats.parquet

    python3 stream_meta_results.py --local \
        --meta-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v1/meta_results \
        --out=full_sumstats.parquet
"""

import io
import os
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

import click
import pandas as pd

# Columns kept from each per-gene meta TSV. Everything else (motif strings, per-cohort allele
# frequency JSON, the random-effects betas/CIs) is dropped here -- it is recoverable from the source
# and would multiply the output size for no gain in the analyses this file feeds. The per-cohort
# betas/SEs (_1 = BioHEART, _2 = TOB) are kept for cross-cohort replication and heterogeneity.
KEEP = [
    'federate_id', 'pair_class', 'cell_type', 'chr', 'gene',
    'bioheart_vid', 'bioheart_pos', 'bioheart_svtype',
    'tob_vid', 'tob_pos', 'tob_svtype',
    'n_samples_tested_1', 'n_samples_tested_2',
    'coeff_1', 'se_1', 'pval_1', 'alt_af_1',
    'coeff_2', 'se_2', 'pval_2', 'alt_af_2',
    'coeff_meta_fixed', 'se_meta_fixed', 'pval_meta_fixed', 'pval_meta_random',
    'pval_q_meta', 'q_meta', 'i2_meta', 'direction_concordant', 'abs_alt_af_diff',
]


def _schema():
    """Output schema. Built lazily so the submitting driver does not need pyarrow, only the jobs."""
    import pyarrow as pa  # noqa: PLC0415

    # Dictionary-encoded strings: these repeat heavily across rows, so dictionary encoding is what
    # keeps the output in the tens of MB rather than hundreds.
    def cat():
        return pa.dictionary(pa.int32(), pa.string())

    return pa.schema([
        ('federate_id', cat()),
        ('pair_class', cat()),
        ('cell_type', cat()),
        ('chr', cat()),
        ('gene_name', cat()),
        ('bioheart_vid', cat()),
        ('bioheart_pos', pa.int32()),
        ('bioheart_svtype', cat()),
        ('tob_vid', cat()),
        ('tob_pos', pa.int32()),
        ('tob_svtype', cat()),
        ('n_samples_tested_1', pa.int16()),
        ('n_samples_tested_2', pa.int16()),
        ('coeff_1', pa.float32()),
        ('se_1', pa.float32()),
        ('pval_1', pa.float64()),
        ('alt_af_1', pa.float32()),
        ('coeff_2', pa.float32()),
        ('se_2', pa.float32()),
        ('pval_2', pa.float64()),
        ('alt_af_2', pa.float32()),
        ('coeff_meta_fixed', pa.float32()),
        ('se_meta_fixed', pa.float32()),
        # p-values stay float64: fixed-effect p reaches 1e-300 and float32 would flush it to zero
        ('pval_meta_fixed', pa.float64()),
        ('pval_meta_random', pa.float64()),
        ('pval_q_meta', pa.float32()),
        ('q_meta', pa.float32()),
        ('i2_meta', pa.float32()),
        ('direction_concordant', pa.bool_()),
        ('abs_alt_af_diff', pa.float32()),
    ])


def _split(gs_path):
    """'gs://bucket/a/b' -> ('bucket', 'a/b')."""
    bucket, _, name = gs_path.replace('gs://', '').rstrip('/').partition('/')
    return bucket, name


def _exists(client, path):
    if not path.startswith('gs://'):
        return os.path.exists(path)
    bucket, name = _split(path)
    return client.bucket(bucket).blob(name).exists()


def _put(client, local, path):
    """Copy a finished local file to path (gs:// or local), so a partial file never appears there."""
    if path.startswith('gs://'):
        bucket, name = _split(path)
        client.bucket(bucket).blob(name).upload_from_filename(local)
        os.remove(local)
    else:
        os.replace(local, path)


def _read_one(blob):
    """Read one gene's meta TSV from GCS, projected to KEEP.

    Returns None if empty, and the exception if unreadable, so the caller can count failures
    instead of silently losing genes.
    """
    try:
        df = pd.read_csv(io.BytesIO(blob.download_as_bytes()), sep='\t')
    except Exception as exc:  # noqa: BLE001 - a missing or truncated gene file must not kill the sweep
        return exc
    if df.empty:
        return None
    return df[[c for c in KEEP if c in df.columns]]


def _list_prefixes(client, bucket_name, prefix):
    """Immediate 'subdirectories' of prefix (one level, like `gsutil ls`)."""
    it = client.list_blobs(bucket_name, prefix=prefix, delimiter='/')
    for _ in it.pages:  # prefixes are only populated once the pages have been consumed
        pass
    return sorted(it.prefixes)


def _to_arrow(df):
    """Coerce a concatenated cell-type frame to the output schema."""
    import pyarrow as pa  # noqa: PLC0415

    schema = _schema()
    df = df.rename(columns={'gene': 'gene_name'})
    for col in ('n_samples_tested_1', 'n_samples_tested_2', 'bioheart_pos', 'tob_pos'):
        df[col] = pd.to_numeric(df.get(col), errors='coerce').fillna(0)
    if 'direction_concordant' in df:
        df['direction_concordant'] = df['direction_concordant'].astype('boolean').fillna(False)
    else:
        df['direction_concordant'] = False
    for field in schema:
        if field.name not in df.columns:
            df[field.name] = None
    return pa.Table.from_pandas(df[[f.name for f in schema]], schema=schema, preserve_index=False)


def stream_cell_type(meta_dir, cell_type, part_path, threads, compression):
    """One cell type's per-gene meta TSVs -> one Parquet part file (skipped if it already exists)."""
    import pyarrow.parquet as pq  # noqa: PLC0415
    from google.cloud import storage  # noqa: PLC0415

    client = storage.Client()
    if _exists(client, part_path):
        print(f'{cell_type}: {part_path} already exists, skipped')
        return

    bucket_name, prefix = _split(meta_dir)
    t = time.time()
    files = [b for b in client.list_blobs(bucket_name, prefix=f'{prefix}/{cell_type}/') if b.name.endswith('.tsv')]
    with ThreadPoolExecutor(threads) as ex:
        results = list(ex.map(_read_one, files))

    failed = [(b.name, r) for b, r in zip(files, results) if isinstance(r, Exception)]
    for name, exc in failed[:5]:
        print(f'  !! unreadable: gs://{bucket_name}/{name}: {exc!r}', flush=True)
    frames = [r for r in results if isinstance(r, pd.DataFrame)]
    if not frames:
        raise RuntimeError(f'{cell_type}: no readable gene files under {meta_dir}/{cell_type}')

    table = _to_arrow(pd.concat(frames, ignore_index=True))
    local = tempfile.NamedTemporaryFile(suffix='.parquet', delete=False).name
    pq.write_table(table, local, compression=compression, row_group_size=100_000)
    _put(client, local, part_path)

    dt = time.time() - t
    print(f'{cell_type}: {len(files):,} files -> {table.num_rows:,} rows in {dt:.0f}s '
          f'({len(files) / max(dt, 1e-9):.0f} files/s) -> {part_path}', flush=True)
    if failed:
        print(f'!! {cell_type}: {len(failed):,} gene files could not be read and are missing', file=sys.stderr)


def combine_parts(part_paths, out, compression):
    """Concatenate the per-cell-type parts into one Parquet, cell type by cell type."""
    import pyarrow.parquet as pq  # noqa: PLC0415
    from google.cloud import storage  # noqa: PLC0415

    client = storage.Client()
    missing = [p for p in part_paths if not _exists(client, p)]
    if missing:
        raise RuntimeError(f'{len(missing)} part file(s) missing, not combining: {missing}')

    schema = _schema()
    local_out = tempfile.NamedTemporaryFile(suffix='.parquet', delete=False).name
    total = 0
    with pq.ParquetWriter(local_out, schema, compression=compression) as writer:
        for part in part_paths:
            if part.startswith('gs://'):
                bucket, name = _split(part)
                with tempfile.NamedTemporaryFile(suffix='.parquet') as tmp:
                    client.bucket(bucket).blob(name).download_to_filename(tmp.name)
                    table = pq.read_table(tmp.name, schema=schema)
            else:
                table = pq.read_table(part, schema=schema)
            # one row group per <=100k rows of a single cell type keeps cell_type filters selective
            writer.write_table(table, row_group_size=100_000)
            total += table.num_rows
            print(f'  {part}: {table.num_rows:,} rows', flush=True)
    _put(client, local_out, out)
    print(f'\n{total:,} rows from {len(part_paths)} cell types -> {out}')
    print('All columns are summary statistics, identifiers, allele frequencies or sample counts.')


@click.command()
@click.option('--meta-dir', required=True, help='gs:// path to meta_results')
@click.option('--out', default='full_sumstats.parquet', help='combined Parquet, local path or gs://')
@click.option('--parts-dir', default=None, help='where per-cell-type parts go; default: {out minus .parquet}_parts')
@click.option('--cell-types', default='', help='comma-separated subset; default is every cell type found')
@click.option('--threads', default=64, type=int, help='concurrent GCS reads per cell type')
@click.option('--compression', default='zstd')
@click.option('--local', is_flag=True, help='run in-process, one cell type at a time, instead of a Hail Batch')
@click.option('--image', default=None, help='image for the batch jobs (needs pandas, pyarrow, google-cloud-storage); default: the driver image')
@click.option('--job-cpu', default=2, type=float, help='CPUs per cell-type job')
def main(meta_dir, out, parts_dir, cell_types, threads, compression, local, image, job_cpu):
    """Stream meta_results from GCS into one compact Parquet, without persisting the TSVs."""
    from google.cloud import storage  # noqa: PLC0415

    client = storage.Client()
    bucket_name, prefix = _split(meta_dir)
    found = [p.rstrip('/').rsplit('/', 1)[-1] for p in _list_prefixes(client, bucket_name, f'{prefix}/')]
    wanted = [c.strip() for c in cell_types.split(',') if c.strip()] or found
    unknown = [c for c in wanted if c not in found]
    if unknown:
        raise click.UsageError(f'cell type(s) not present under {meta_dir}: {unknown}')

    parts_dir = (parts_dir or out.removesuffix('.parquet') + '_parts').rstrip('/')
    part_paths = [f'{parts_dir}/{ct}.parquet' for ct in wanted]
    print(f'{len(wanted)} cell types; parts -> {parts_dir}/; combined -> {out}')

    if local:
        if not parts_dir.startswith('gs://'):
            os.makedirs(parts_dir, exist_ok=True)
        t0 = time.time()
        for ct, part in zip(wanted, part_paths):
            stream_cell_type(meta_dir, ct, part, threads, compression)
        combine_parts(part_paths, out, compression)
        print(f'done in {(time.time() - t0) / 60:.1f} min')
        return

    from cpg_utils.hail_batch import get_batch  # noqa: PLC0415

    b = get_batch(name='stream SV meta results to parquet')
    jobs = []
    for ct, part in zip(wanted, part_paths):
        j = b.new_python_job(name=f'stream_{ct}')
        if image:
            j.image(image)
        j.cpu(job_cpu)
        j.call(stream_cell_type, meta_dir, ct, part, threads, compression)
        jobs.append(j)

    combine = b.new_python_job(name='combine_parts')
    if image:
        combine.image(image)
    combine.depends_on(*jobs)
    combine.call(combine_parts, part_paths, out, compression)

    b.run(wait=False)


if __name__ == '__main__':
    main()
