#!/usr/bin/env python3
"""
Hail Batch runner for sv_cell_type_specificity.py (flagship cell-type specificity scenarios, v2 SV meta-analysis).

Two jobs:
    gather    Python job (default image: pandas, pyarrow, google-cloud-storage). Writes three small TSVs to
              {out}/inputs/: the gene-level q-values of every cell type, the expressed-gene lists (either cohort), and
              every SV test of every eGene's gene (all cell types) from the v2 full_sumstats.parquet.
    classify  trtools image (trtools/cyvcf2 for SV dosages). Gets those TSVs and the 44 per-cohort SV VCFs as Batch
              inputs and runs sv_cell_type_specificity.py on local files (no GCS access inside the job).

Outputs in {out}/: esv_cell_type_specificity.tsv, sv_pair_ld.tsv, rescue_tests.tsv (see sv_cell_type_specificity.py).

analysis-runner --dataset tenk10k --access-level test --description "SV cell-type specificity (flagship scenarios, v2)" \
    --output-dir tenk10k-sv/sv/cell_type_spec/bioheart_n968_and_tob_n935/7pc/v2 \
    python3 cell_type_spec_runner.py
"""
import os

import click
from cpg_utils.config import get_config

HERE = os.path.dirname(os.path.abspath(__file__))
V2 = 'gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v2_F_filtered'
GENE_LIST_ROOT = 'gs://cpg-tenk10k-test/sv/input_files/{cohort}/7pc/scRNA_gene_lists/1_min_pct_cells_expressed'
SV_VCF = ('gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/'
          '{cohort}_common_chr_specific/hail_filtered_{chrom}.vcf.bgz')
COHORT_DIRS = {'bioheart': 'bioheart_n968', 'tob': 'tob_n935'}
SS_COLS = ['gene_name', 'cell_type', 'chr', 'federate_id', 'pair_class', 'bioheart_vid', 'bioheart_pos', 'tob_vid',
           'coeff_meta_fixed', 'se_meta_fixed']


def gather(qval_dir, gene_list_roots, sumstats, out_q, out_expr, out_ss):
    """Self-contained (runs as a Batch Python job): q-values, expressed genes and the eGenes' SV tests -> 3 TSVs."""
    import io  # noqa: PLC0415
    import json  # noqa: PLC0415
    import tempfile  # noqa: PLC0415
    from concurrent.futures import ThreadPoolExecutor  # noqa: PLC0415
    import pandas as pd  # noqa: PLC0415
    from google.cloud import storage  # noqa: PLC0415

    client = storage.Client()

    def split(p):
        bucket, _, name = p.removeprefix('gs://').partition('/')
        return bucket, name

    def blobs(prefix, suffix):
        bucket, name = split(prefix)
        return [(bucket, b.name) for b in client.list_blobs(bucket, prefix=name.rstrip('/') + '/') if b.name.endswith(suffix)]

    def get(bn):
        return client.bucket(bn[0]).blob(bn[1]).download_as_bytes()

    def put(df, path):
        if not path.startswith('gs://'):   # local test
            df.to_csv(path, sep='\t', index=False, compression='gzip')
            return
        bucket, name = split(path)
        with tempfile.NamedTemporaryFile(suffix='.tsv.gz') as fh:
            df.to_csv(fh.name, sep='\t', index=False, compression='gzip')
            client.bucket(bucket).blob(name).upload_from_filename(fh.name)

    def read_q(bn):
        d = pd.read_csv(io.BytesIO(get(bn)), sep='\t', usecols=['gene_name', 'qval'])
        return d.assign(cell_type=bn[1].rsplit('/', 1)[1].removesuffix('_qval.tsv'))

    with ThreadPoolExecutor(16) as ex:
        q = pd.concat(ex.map(read_q, blobs(qval_dir, '_qval.tsv')), ignore_index=True)
    put(q, out_q)
    print(f'q-values: {len(q):,} gene x cell type rows; {int((q.qval < 0.05).sum()):,} at q < 0.05')

    def read_list(bn):
        ct = bn[1].rstrip('/').split('/')[-2]
        return [(g, ct) for g in json.loads(get(bn))]

    pairs = set()
    with ThreadPoolExecutor(32) as ex:
        for root in gene_list_roots:
            for chunk in ex.map(read_list, blobs(root, '_gene_list.json')):
                pairs.update(chunk)
    put(pd.DataFrame(sorted(pairs), columns=['gene_name', 'cell_type']), out_expr)
    print(f'expressed gene x cell type pairs (either cohort): {len(pairs):,}')

    with tempfile.NamedTemporaryFile(suffix='.parquet') as fh:
        client.bucket(split(sumstats)[0]).blob(split(sumstats)[1]).download_to_filename(fh.name)
        ss = pd.read_parquet(fh.name, columns=SS_COLS)
    for k in ss.select_dtypes('category').columns:
        ss[k] = ss[k].astype(str)
    egenes = set(q.loc[q.qval < 0.05, 'gene_name'])
    ss = ss[ss.gene_name.isin(egenes)]
    put(ss, out_ss)
    print(f'SV tests of the {len(egenes):,} eGenes, all cell types: {len(ss):,} rows')


@click.command()
@click.option('--qval-dir', default=f'{V2}/fdr_qvals/using_acat')
@click.option('--sumstats', default=f'{V2}/full_sumstats.parquet')
@click.option('--rescue-effect', default='abs', type=click.Choice(['abs', 'signed']))
@click.option('--memory', default='highmem')
def main(qval_dir, sumstats, rescue_effect, memory):
    from cpg_utils.hail_batch import get_batch, output_path  # noqa: PLC0415

    b = get_batch(name='SV cell-type specificity (flagship scenarios)')
    inputs = {k: output_path(f'inputs/{k}.tsv.gz', 'analysis') for k in ('qvals', 'expressed', 'sumstats')}

    g = b.new_python_job(name='gather inputs')
    g.memory('highmem')
    g.call(gather, qval_dir, [GENE_LIST_ROOT.format(cohort=c) for c in COHORT_DIRS.values()], sumstats,
           inputs['qvals'], inputs['expressed'], inputs['sumstats'])

    srcs = {n: open(os.path.join(HERE, *p)).read() for n, p in {
        'sv_cell_type_specificity.py': ('sv_cell_type_specificity.py',),
        'finemap_prep.py': ('..', 'finemap_mcv', 'finemap_prep.py'),
        'fast_associatr.py': ('..', 'snp_tr_1mb', 'fast_associatr.py')}.items()}
    local = {k: b.read_input(v) for k, v in inputs.items()}
    vcfs = {(c, f'chr{i}'): b.read_input(SV_VCF.format(cohort=c, chrom=f'chr{i}'))
            for c in COHORT_DIRS for i in range(1, 23)}

    j = b.new_job(name='classify scenarios')
    j.depends_on(g)
    j.image(get_config()['images']['trtools'])
    j.memory(memory)
    j.storage('10G')
    j.declare_resource_group(out={'esv': '{root}.esv.tsv', 'ld': '{root}.ld.tsv', 'rescue': '{root}.rescue.tsv'})
    j.command(
        'set -euo pipefail\nmkdir -p work/vcfs && cd work\n'
        + ''.join(f"cat > {n} <<'__SRC_EOF__'\n{s}\n__SRC_EOF__\n" for n, s in srcs.items())
        + ''.join(f'ln -s {v} vcfs/{c}_{chrom}.vcf.bgz\n' for (c, chrom), v in vcfs.items())
        + f"python3 sv_cell_type_specificity.py --sumstats {local['sumstats']} --qvals {local['qvals']} "
          f"--expressed {local['expressed']} --vcf-cache vcfs --out-dir out --rescue-effect {rescue_effect} "
          '2>&1 | { grep -v -i "warn" || true; }\n'
        + f"mv out/esv_cell_type_specificity.tsv {j.out['esv']}\n"
        + f"mv out/sv_pair_ld.tsv {j.out['ld']}\n"
        + f"mv out/rescue_tests.tsv {j.out['rescue']}\n",
    )
    b.write_output(j.out['esv'], output_path('esv_cell_type_specificity.tsv', 'analysis'))
    b.write_output(j.out['ld'], output_path('sv_pair_ld.tsv', 'analysis'))
    b.write_output(j.out['rescue'], output_path('rescue_tests.tsv', 'analysis'))
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
