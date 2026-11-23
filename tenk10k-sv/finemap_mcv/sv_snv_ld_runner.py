#!/usr/bin/env python3
# pylint: disable=import-error,wrong-import-position,too-many-arguments
"""
SV-SNV LD for every lead sc-eSV: for each lead SV pair (one per eGene x cell type at gene-level q < 0.05 in the
v2 F-filtered SV meta-analysis), r^2 with every small variant within +-1 Mb, pooled over BioHEART + TOB exactly as
the pooled fine-mapping model (see sv_snv_ld.py). Used to quantify how resolvable each SV is from tagging SNVs
(finemap_results.ipynb, section 16).

Jobs: one per chromosome (trtools image: cyvcf2, numpy, pandas). Each reads the jointly called common-variant SNV
VCF and both cohorts' SV VCFs end to end; the largest chromosomes hold ~1M variants x 1,903 samples (~2 GB int8), so
the default is a 'highmem' 1-CPU job (~6.5 GiB).

Output: {out}/sv_snv_ld/{chrom}.tsv (analysis bucket). Chromosomes with an existing output are skipped, so a failed
batch can be resubmitted with the same command.

analysis-runner --dataset tenk10k --access-level test --description "SV-SNV LD for lead sc-eSVs (v2)" \
    --output-dir tenk10k-sv/sv/finemap_mcv/bioheart_n968_and_tob_n935/7pc/v2 \
    python3 sv_snv_ld_runner.py

Test on one small chromosome first: add --chromosomes=chr21
"""
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.config import get_config

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'snp_tr_1mb'))
from snp_tr_cis_runner import SNP_VCF  # noqa: E402

SV_VCF = ('gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/'
          '{cohort}_common_chr_specific/hail_filtered_{chrom}.vcf.bgz')
QVAL_DIR = ('gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/'
            'v2_F_filtered/fdr_qvals/using_acat')
CHROMS = ','.join(f'chr{i}' for i in range(1, 23))


def lead_pairs(qval_dir):
    """Distinct lead SV pairs of all eGene x cell types at q < 0.05. Tied leads are stored as comma-separated,
    position-aligned lists in bioheart_vid / tob_vid, so they are split in parallel."""
    files = [f for f in to_path(qval_dir).glob('*_qval.tsv')]

    def read(f):
        with f.open() as h:
            return pd.read_csv(h, sep='\t', usecols=['gene_name', 'qval', 'bioheart_vid', 'tob_vid'])

    with ThreadPoolExecutor(16) as ex:
        q = pd.concat(ex.map(read, files), ignore_index=True)
    q = q[q.qval < 0.05]
    rows, skipped = [], 0
    for b, t in zip(q.bioheart_vid.astype(str), q.tob_vid.astype(str)):
        bs, ts = b.split(','), t.split(',')
        if len(bs) != len(ts):
            skipped += 1
            continue
        rows.extend(zip(bs, ts))
    pairs = pd.DataFrame(rows, columns=['bioheart_vid', 'tob_vid']).drop_duplicates()
    pairs['chrom'] = pairs.bioheart_vid.str.extract(r'all_batches\.(chr[0-9XY]+)\.', expand=False)
    print(f'{len(q):,} eGene x cell types -> {len(pairs):,} distinct lead SV pairs '
          f'({skipped} rows with mismatched tie lists skipped; {pairs.chrom.isna().sum()} without a chromosome)')
    return pairs.dropna(subset=['chrom'])


@click.command()
@click.option('--qval-dir', default=QVAL_DIR)
@click.option('--chromosomes', default=CHROMS)
@click.option('--window', default=1_000_000, type=int, help='bp either side of the SV')
@click.option('--memory', default='highmem', help="Batch memory ('standard', 'highmem' or e.g. '10Gi')")
@click.option('--cpu', default=1, type=int)
@click.option('--storage', default='10G')
@click.option('--spot/--no-spot', default=True)
def main(qval_dir, chromosomes, window, memory, cpu, storage, spot):
    from cpg_utils.hail_batch import get_batch, output_path

    pairs = lead_pairs(qval_dir)
    src = open(os.path.join(HERE, 'sv_snv_ld.py')).read()
    b = get_batch(name='SV-SNV LD for lead sc-eSVs')
    n = 0
    for chrom in chromosomes.split(','):
        out = output_path(f'sv_snv_ld/{chrom}.tsv', 'analysis')
        if to_path(out).exists():
            print(f'{chrom}: output exists, skipped')
            continue
        leads = pairs[pairs.chrom == chrom][['bioheart_vid', 'tob_vid']]
        if leads.empty:
            continue
        snp = b.read_input(SNP_VCF.format(chrom=chrom))
        bh = b.read_input(SV_VCF.format(cohort='bioheart', chrom=chrom))
        tob = b.read_input(SV_VCF.format(cohort='tob', chrom=chrom))
        j = b.new_job(name=f'SV-SNV LD {chrom}')
        j.image(get_config()['images']['trtools'])
        j.cpu(cpu)
        j.memory(memory)
        j.storage(storage)
        if hasattr(j, 'spot'):
            j.spot(spot)
        j.command(
            'set -euo pipefail\nmkdir -p work && cd work\n'
            f"cat > sv_snv_ld.py <<'__LD_EOF__'\n{src}\n__LD_EOF__\n"
            f"cat > leads.tsv <<'__LEADS_EOF__'\n{leads.to_csv(sep=chr(9), index=False)}__LEADS_EOF__\n"
            f'python3 sv_snv_ld.py {snp} {bh} {tob} leads.tsv {j.ofile} {window} 2>&1 | '
            '{ grep -v -i "warn\\|bcf_hdr\\|vcf_parse" || true; }\n'
            f'test -s {j.ofile}\n',
        )
        b.write_output(j.ofile, out)
        n += 1
        print(f'{chrom}: {len(leads):,} lead SV pairs')
    print(f'{n} jobs')
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
