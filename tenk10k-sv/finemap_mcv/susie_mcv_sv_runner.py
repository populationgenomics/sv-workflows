#!/usr/bin/env python3
# pylint: disable=too-many-locals,too-many-arguments,import-error,wrong-import-position
"""
SuSiE multiple-causal-variant (MCV) fine-mapping of SNPs + TRs + SVs: per cohort, pooled, or both.

Scope: the 30,799 eGene x cell type pairs (8,111 eGenes) with gene-level q < 0.05 in the F-filtered v2 SV
meta-analysis (7pc/v2_F_filtered: fixed FederateSV pairing, recovered pairs included), in each pair's own cell type
only. --qval-dir / --meta-set select another meta-analysis (they must come from the same run).

Two designs, run from one shared prep (--mode both, the default)
---------------------------------------------------------------
SNPs and TRs were called jointly across BioHEART and TOB; SVs were called per cohort and paired afterwards by
FederateSV (V2_META_SET).

  per_cohort  each cohort on its own (n ~ 950): every SV in that cohort's callset is a candidate, LD is exactly
              in-sample, and the two cohorts give a replication test for every credible set (combine_cohorts.py).
  pooled      BioHEART + TOB in one model (n = 1,903): candidates are the variants the meta-analysis could pool --
              SNPs/TRs, and federated SV pairs (concordant_* and dup_cnv, never del_cnv) -- each passing
              associaTR's filters in both cohorts. A federated SV is one column holding each cohort's own dosage.
              Imputation and covariate residualisation happen within each cohort before stacking. Every
              reported eSV is a federated pair, so the pooled model covers all of them. See susie_mcv.R.

Inputs (identical to the SV eQTL scan; v2 changed only the pairing, so the per-cohort 'v1' input folders are current):
    cis window   {SV_INPUTS}/bioheart_n968/7pc/cis_window_files/v1/{ct}/{chrom}/{gene}_1000000bp.bed
    pheno_cov    {SV_INPUTS}/{cohort_dir}/7pc/pheno_cov_numpy/v1/{ct}/{chrom}/{gene}_pheno_cov.npy
    SNP / TR     the jointly-called SNP and TR VCFs used by snp_tr_cis_runner.py
    SV           the per-cohort EH-style SV VCFs used by associatr_runner_sv.py

Jobs: per (chromosome, chunk of --genes-per-job genes)
    prep   trtools image: finemap_prep.py builds dosage matrices for both cohorts -> one tarball
    susie  R image with susieR: susie_mcv.R runs every gene x cell type in the chunk, per cohort and/or pooled
Outputs:
    {out}/variants/{chrom}_{chunk:03d}.tsv.gz   PIP, CS membership, posterior effect (PIP >= 0.01, CS members, all SVs)
    {out}/cs/{chrom}_{chunk:03d}.tsv.gz         one row per credible set
    {out}/runs/{chrom}_{chunk:03d}.tsv.gz       one row per gene x cohort x cell type (cohort = bioheart/tob/pooled)
    {out}/prep_log/{chrom}_{chunk:03d}.tsv      samples and variants per gene x cohort
Existing cs/ outputs are skipped, so a failed batch can be resubmitted.

Pilot (10 genes on chr22, one prep + one SuSiE job):

analysis-runner --dataset tenk10k --access-level test --description "SuSiE MCV SNP+TR+SV pilot (v2)" \
    --output-dir tenk10k-sv/sv/finemap_mcv/bioheart_n968_and_tob_n935/7pc/v2_pilot \
    python3 susie_mcv_sv_runner.py --chromosomes=chr22 --max-genes=10 --genes-per-job=10

Full run (resubmit the same command to fill in failed chunks):

analysis-runner --dataset tenk10k --access-level test --description "SuSiE MCV SNP+TR+SV per cohort and pooled (v2)" \
    --output-dir tenk10k-sv/sv/finemap_mcv/bioheart_n968_and_tob_n935/7pc/v2 \
    python3 susie_mcv_sv_runner.py

Long or preempted chunks: resubmit with --split=5 (keep --genes-per-job=25) to rerun each unfinished chunk as 5-gene
sub-jobs; add --no-spot to stop preemption. Cancel the old batch first, or a chunk could finish in both.

Jobs killed with exit 137 ran out of memory (1 CPU 'standard' is ~3.75 GiB). Resubmit the same command with more
memory: only chunks without a cs/ output are rebuilt, e.g. --susie-memory=highmem --prep-memory=highmem (~6.5 GiB at
1 CPU), or add --susie-cpu=2 / --prep-cpu=2 to double that again.
"""
import json
import os
import sys

import click

from cpg_utils import to_path
from cpg_utils.config import get_config

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'snp_tr_1mb'))
from snp_tr_cis_runner import SNP_VCF, SV_INPUTS, TR_VCF, cis_region, egene_pairs, gene_chrom  # noqa: E402

COHORTS = {'bioheart': 'bioheart_n968', 'tob': 'tob_n935'}
SV_VCF = 'gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/{cohort}_common_chr_specific/hail_filtered_{chrom}.vcf.bgz'
R_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/r-meta:r-fast'   # used by str/fine-mapping/mcv
# v2 meta-analysis and the pairing table it pooled. The table is federate_overlap_out/meta_analysis_set_fixed.csv
# restricted to the F filter (15,109 pairs; recovered rows whose federate_id collided carry a __row<N> suffix).
# Checked against v2_F_filtered/full_sumstats.parquet: all 14,424 tested pairs are in it; the other 685 were never
# tested (SV not in an eQTL VCF, failed filters, or no cis gene), and drop out of the pooled model the same way.
V2_META = ('gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/'
           'v2_F_filtered')
V2_META_SET = 'gs://cpg-tenk10k-test/sv/input_files/meta_analysis_set_v2_fixed_F_filtered.csv'
# For v1_F_filtered, the matching table is the original meta_analysis_set.csv restricted to F-passing pairs
# (13,576 pairs, meta_analysis_set.as_used_v1_F_filtered.csv; not uploaded).


@click.command()
@click.option('--qval-dir', default=f'{V2_META}/fdr_qvals/using_acat')
@click.option('--cohorts', default='bioheart,tob', help='comma-separated; each is fine-mapped on its own')
@click.option('--genes-per-job', default=25, type=int,
              help='genes per chunk; defines the output tags, so keep it fixed when resubmitting (use --split)')
@click.option('--split', default=0, type=int,
              help='run each unfinished chunk as sub-jobs of at most this many genes (0 = whole chunk)')
@click.option('--spot/--no-spot', default=True, help='--no-spot: non-preemptible VMs (dearer, never preempted)')
@click.option('--chromosomes', default=','.join(f'chr{i}' for i in range(1, 23)))
@click.option('--max-genes', default=0, type=int, help='pilot: cap genes per chromosome (0 = all)')
@click.option('--max-parallel-jobs', default=300, type=int)
@click.option('--num-causal', default=10, type=int, help='SuSiE L')
@click.option('--max-iter', default=100, type=int)
@click.option('--mode', default='both', type=click.Choice(['per_cohort', 'pooled', 'both']),
              help='per-cohort models, the pooled BioHEART+TOB model, or both (shared prep)')
@click.option('--meta-set', default=V2_META_SET,
              help='pairing table the meta-analysis in --qval-dir actually pooled (see V2_META_SET)')
@click.option('--r-image', default=None, help=f'image with Rscript + susieR (default: config images.susie, else {R_IMAGE})')
@click.option('--prep-storage', default='15G')
@click.option('--prep-memory', default='standard', help="Batch memory for prep jobs ('standard', 'highmem' or e.g. '8Gi')")
@click.option('--prep-cpu', default=1, type=int)
@click.option('--susie-memory', default='standard')
@click.option('--susie-cpu', default=1, type=int, help='susieR is single-threaded')
def main(qval_dir, cohorts, genes_per_job, split, spot, chromosomes, max_genes, max_parallel_jobs, num_causal, max_iter, mode,
         meta_set, r_image,
         prep_storage, prep_memory, prep_cpu, susie_memory, susie_cpu):
    from cpg_utils.hail_batch import get_batch, output_path

    cohorts = cohorts.split(',')
    assert set(cohorts) <= set(COHORTS), cohorts
    if mode != 'per_cohort':
        assert set(cohorts) == set(COHORTS), '--mode pooled/both needs both cohorts'
    r_image = r_image or get_config().get('images', {}).get('susie') or R_IMAGE
    pairs = egene_pairs(qval_dir)
    pairs['chrom'] = pairs.gene_name.map(gene_chrom(pairs))
    assert pairs.chrom.notna().all(), pairs[pairs.chrom.isna()].head()
    print(f'{len(pairs):,} eGene x cell type pairs, {pairs.gene_name.nunique():,} genes; cohorts {cohorts}')

    srcs = {n: open(os.path.join(HERE, *p)).read() for n, p in {
        'fast_associatr.py': ('..', 'snp_tr_1mb', 'fast_associatr.py'),
        'finemap_prep.py': ('finemap_prep.py',), 'susie_mcv.R': ('susie_mcv.R',)}.items()}

    b = get_batch(name=f'SuSiE MCV SNP+TR+SV ({mode})')
    if mode != 'per_cohort' and not to_path(meta_set).exists():
        raise click.UsageError(f'pairing table {meta_set} not found (see V2_META_SET)')
    meta_set_file = b.read_input(meta_set) if mode != 'per_cohort' else None
    jobs, n_pairs = [], 0

    done_tags = {f.name.removesuffix('.tsv.gz') for f in to_path(output_path('cs', 'analysis')).glob('*.tsv.gz')}

    def work_units(chrom, genes):
        """(tag, genes) still to run. A chunk is genes[k*genes_per_job : (k+1)*genes_per_job], tagged {chrom}_{k:03d}.
        With --split, an unfinished chunk runs as sub-jobs tagged {chrom}_{k:03d}_g{a:02d}-{b:02d} (gene offsets a..b
        within the chunk, inclusive). Offsets covered by existing sub-job outputs are skipped, so --split can be changed
        between resubmissions without duplicating or missing genes."""
        for chunk, start in enumerate(range(0, len(genes), genes_per_job)):
            tag = f'{chrom}_{chunk:03d}'
            block = genes[start:start + genes_per_job]
            if tag in done_tags:
                continue
            covered = set()
            for t in done_tags:
                if t.startswith(f'{tag}_g'):
                    a, b = map(int, t.removeprefix(f'{tag}_g').split('-'))
                    covered.update(range(a, b + 1))
            if not split and not covered:
                yield tag, block
                continue
            todo = [i for i in range(len(block)) if i not in covered]
            runs = []   # contiguous runs of offsets, each cut into pieces of <= split genes
            for i in todo:
                if runs and i == runs[-1][-1] + 1 and len(runs[-1]) < (split or len(block)):
                    runs[-1].append(i)
                else:
                    runs.append([i])
            for r in runs:
                yield f'{tag}_g{r[0]:02d}-{r[-1]:02d}', [block[i] for i in r]

    def throttle(j):
        if hasattr(j, 'spot'):
            j.spot(spot)
        if len(jobs) >= max_parallel_jobs:
            j.depends_on(jobs[-max_parallel_jobs])
        jobs.append(j)

    for chrom in chromosomes.split(','):
        genes = sorted(pairs.loc[pairs.chrom == chrom, 'gene_name'].unique())
        if max_genes:
            genes = genes[:max_genes]
        if not genes:
            continue
        snp = b.read_input_group(base=SNP_VCF.format(chrom=chrom), tbi=SNP_VCF.format(chrom=chrom) + '.tbi')
        tr = b.read_input_group(base=TR_VCF.format(chrom=chrom), tbi=TR_VCF.format(chrom=chrom) + '.tbi')
        # read end to end by finemap_prep.py, so no index is needed (the TOB SV .tbi files have been stale before)
        svs = {c: b.read_input(SV_VCF.format(cohort=c, chrom=chrom)) for c in cohorts}
        for tag, chunk_genes in work_units(chrom, genes):
            manifest = {'snp_vcf': str(snp.base), 'tr_vcf': str(tr.base),
                        'sv_vcf': {c: str(svs[c]) for c in cohorts}, 'genes': []}
            if meta_set_file is not None:
                manifest['meta_set'] = str(meta_set_file)
            for gene in chunk_genes:
                cts = sorted(pairs.loc[pairs.gene_name == gene, 'cell_type'])
                n_pairs += len(cts)
                manifest['genes'].append({
                    'gene': gene, 'region': cis_region(cts[0], chrom, gene),
                    'pheno': {ct: {c: str(b.read_input(
                        f'{SV_INPUTS}/{COHORTS[c]}/7pc/pheno_cov_numpy/v1/{ct}/{chrom}/{gene}_pheno_cov.npy'))
                        for c in cohorts} for ct in cts}})

            prep = b.new_job(name=f'finemap prep {tag}')
            prep.image(get_config()['images']['trtools'])
            prep.storage(prep_storage)
            prep.cpu(prep_cpu)
            prep.memory(prep_memory)
            prep.command(
                'set -euo pipefail\nmkdir -p work && cd work\n'
                f"cat > fast_associatr.py <<'__FAST_EOF__'\n{srcs['fast_associatr.py']}\n__FAST_EOF__\n"
                f"cat > finemap_prep.py <<'__PREP_EOF__'\n{srcs['finemap_prep.py']}\n__PREP_EOF__\n"
                f"cat > manifest.json <<'__MAN_EOF__'\n{json.dumps(manifest)}\n__MAN_EOF__\n"
                'mkdir -p prep\n'
                'python3 finemap_prep.py manifest.json prep 2>&1 | { grep -v -i "warn" || true; }\n'
                f'cp prep/prep_log.tsv {prep.log}\n'
                f'tar -cf {prep.tar} prep\n',
            )
            b.write_output(prep.log, output_path(f'prep_log/{tag}.tsv', 'analysis'))
            throttle(prep)

            susie = b.new_job(name=f'susie {tag}')
            susie.image(r_image)
            susie.cpu(susie_cpu)
            susie.memory(susie_memory)
            susie.storage(prep_storage)
            susie.declare_resource_group(out={'variants': '{root}.variants.tsv.gz', 'cs': '{root}.cs.tsv.gz',
                                              'runs': '{root}.runs.tsv.gz'})
            susie.command(
                'set -euo pipefail\nmkdir -p work && cd work\n'
                f"cat > susie_mcv.R <<'__R_EOF__'\n{srcs['susie_mcv.R']}\n__R_EOF__\n"
                f'tar -xf {prep.tar}\n'
                f'Rscript susie_mcv.R prep out {num_causal} {max_iter} 0.01 {mode}\n'
                f"mv out.variants.tsv.gz {susie.out['variants']}\n"
                f"mv out.cs.tsv.gz {susie.out['cs']}\n"
                f"mv out.runs.tsv.gz {susie.out['runs']}\n",
            )
            for kind in ('variants', 'runs', 'cs'):   # cs last: its existence marks the chunk as done
                b.write_output(susie.out[kind], output_path(f'{kind}/{tag}.tsv.gz', 'analysis'))
            throttle(susie)
    print(f'{len(jobs)} jobs, {n_pairs:,} gene x cell type pairs x {len(cohorts)} cohort(s)')
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
