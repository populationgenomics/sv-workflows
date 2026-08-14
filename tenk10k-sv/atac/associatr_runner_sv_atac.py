#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments,too-many-locals
"""
Batched SV associaTR runner for ATAC peaks.

associaTR takes one phenotype per invocation, so N peaks always means N invocations. This
script batches peaks_per_job of them into a single pair of Hail Batch jobs to amortise the
per-job overhead, which dominates cost at ~400k peaks per cell type.

Differences vs. str/atac/associatr_runner.py (one job per peak):
- peaks_per_job peaks share one job pair, so the per-chromosome VCF is localised once per
  chunk instead of once per peak.
- the resume check is a single listing of the results directory per cell type, rather than
  one GCS stat per peak.
- results are concatenated in-job, with the peak id added as the first column, so each
  chunk writes one TSV instead of one TSV per peak.

Cis window (as in tenk10k-sv/associatr_runner_sv.py on the associatr_SV-runner branch):
an SV is tested against a peak if its interval [POS, INFO/END] overlaps the peak's cis
window [peak_start - cis_window_size, peak_end + cis_window_size] by >=1bp, rather than
only SVs whose POS falls in the window. associaTR cannot express this, so bcftools
pre-filters per peak and associaTR is then pointed at the whole surviving VCF.

Unlike the branch version the filter also keeps records with no INFO/END (BND, and often
INS), which an 'INFO/END>=start' test alone would silently drop.

A peak whose associaTR call fails will fail its whole chunk (no output is written, so the
chunk is retried on the next run). Failing peaks are all listed in the job log rather than
stopping at the first.

max_peaks caps the peaks run per cell type for smoke tests, taken as a position-ordered
prefix so the chunks it produces match those of an uncapped run; set it to 0 to run all.

 analysis-runner --dataset "bioheart" --config associatr_runner_sv_atac.toml \
    --description "run associatr on SVs for atac peaks" \
    --access-level "test" \
    --output-dir "str/associatr-atac/sv/1mb" \
    python3 associatr_runner_sv_atac.py

"""

from collections import defaultdict

import hail as hl
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config, output_path
from cpg_utils.hail_batch import get_batch, init_batch, reset_batch


def parse_peak(loc):
    """'chr1:65002192-65003572' -> ('chr1', 65002192, 65003572)"""
    chrom, _, span = loc.partition(':')
    start, _, end = span.partition('-')
    return chrom, int(start), int(end)


def main():
    """
    Run associaTR on SVs for ATAC peaks, batching peaks_per_job peaks per job pair
    """
    init_batch()  # for hl.get_reference below
    config = get_config()['associatr_sv_atac']
    window = config['cis_window_size']
    peaks_per_job = config['peaks_per_job']

    for cell_type in config['cell_types'].split(','):
        pheno_dir = f'{config["pheno_cov_numpy_dir"]}/{cell_type}/pheno_cov_numpy'

        # one listing rather than a stat per peak - at 400k peaks the latter dominates driver runtime
        peak_names = [p.name.removesuffix('_pheno_cov.npy') for p in to_path(pheno_dir).glob('*_pheno_cov.npy')]
        results_dir = output_path(f'results/{cell_type}', 'analysis')
        done = {p.name.removesuffix('.tsv') for p in to_path(results_dir).glob('*.tsv')}
        print(f'{cell_type}: {len(peak_names)} peaks with phenotype files, {len(done)} chunks already done')

        # sorted by position so each chunk covers a contiguous span, which keeps the union subset small.
        # max_peaks then takes a positional prefix, so a capped test run produces exactly the same chunks
        # (and therefore the same output names) as the corresponding chunks of an uncapped run.
        peaks = sorted((*parse_peak(name), name) for name in peak_names)
        max_peaks = config.get('max_peaks', 0)
        if max_peaks:
            peaks = peaks[:max_peaks]
            print(f'{cell_type}: capped to {len(peaks)} peaks by max_peaks')

        by_chrom = defaultdict(list)
        for chrom, start, end, name in peaks:
            by_chrom[chrom].append((start, end, name))

        for chrom, entries in sorted(by_chrom.items()):
            reset_batch()
            b = get_batch(name=f'associaTR SV-ATAC {cell_type} {chrom}')
            chrom_len = hl.get_reference('GRCh38').lengths[chrom]
            vcf_path = f'{config["vcf_file_dir"]}/hail_filtered_{chrom}.vcf.bgz'
            variant_vcf = b.read_input_group(base=vcf_path, tbi=f'{vcf_path}.tbi')

            _dependent_jobs: list[hb.batch.job.Job] = []

            def manage_concurrency_for_job(job: hb.batch.job.Job, jobs=_dependent_jobs):
                if len(jobs) >= config['max_parallel_jobs']:
                    job.depends_on(jobs[-config['max_parallel_jobs']])
                jobs.append(job)

            submitted = 0
            for i in range(0, len(entries), peaks_per_job):
                chunk = entries[i : i + peaks_per_job]
                out_name = f'{chrom}_chunk{i // peaks_per_job:05d}_{window}bp'
                if out_name in done:
                    continue

                cis_windows = [(name, max(0, start - window), end + window) for start, end, name in chunk]
                union_start = min(cis_start for _, cis_start, _ in cis_windows)
                union_end = max(cis_end for _, _, cis_end in cis_windows)
                windows_tsv = '\n'.join(f'{name}\t{cis_start}\t{cis_end}' for name, cis_start, cis_end in cis_windows)

                # Job 1: per-peak bcftools overlap filter. bcftools lives in its own image because
                # the trtools image does not ship bcftools/htslib, so the filtered VCFs are handed
                # to job 2 as a tarball.
                filter_job = b.new_job(name=f'Filter SV VCF: {cell_type} {out_name}')
                filter_job.image(get_config()['images']['bcftools'])
                filter_job.cpu(config['filter_job_cpu'])
                filter_job.storage(config['filter_job_storage'])
                filter_job.command(
                    f"""set -euo pipefail
mkdir -p filtered
cat > windows.tsv <<'WINDOWS'
{windows_tsv}
WINDOWS
# subset once to the chunk's union window so the per-peak filters below scan a small file
bcftools view -r {chrom}:1-{union_end} -i 'INFO/END>={union_start} || POS>={union_start}' \
  {variant_vcf.base} -Oz -o chunk.vcf.bgz
tabix -f -p vcf chunk.vcf.bgz
while IFS=$'\\t' read -r peak cis_start cis_end; do
  # [POS, INFO/END] must overlap [cis_start, cis_end]; the POS clause retains records with no END
  bcftools view -r {chrom}:1-${{cis_end}} -i "INFO/END>=${{cis_start}} || POS>=${{cis_start}}" \
    chunk.vcf.bgz -Oz -o filtered/${{peak}}.vcf.bgz
  tabix -f -p vcf filtered/${{peak}}.vcf.bgz
done < windows.tsv
tar -cf {filter_job.ofile} -C filtered .
""",
                )

                # Job 2: associaTR once per peak against its pre-filtered VCF, results concatenated
                # with the peak id prepended (associaTR output does not identify the phenotype).
                assoc_job = b.new_job(name=f'associaTR: {cell_type} {out_name}')
                assoc_job.image(get_config()['images'][config['trtools_image']])
                assoc_job.cpu(config['job_cpu'])
                assoc_job.storage(config['job_storage'])
                assoc_job.declare_resource_group(association_results={'tsv': '{root}.tsv'})

                pheno_tsv = '\n'.join(
                    f'{name}\t{b.read_input(f"{pheno_dir}/{name}_pheno_cov.npy")}' for name, _, _ in cis_windows
                )
                assoc_job.command(
                    f"""set -euo pipefail
mkdir -p filtered results
tar -xf {filter_job.ofile} -C filtered
cat > pheno.tsv <<PHENO
{pheno_tsv}
PHENO
fail=0
while IFS=$'\\t' read -r peak pheno; do
  if ! associaTR results/${{peak}}.tsv filtered/${{peak}}.vcf.bgz ${{peak}} ${{pheno}} \
    --region={chrom}:1-{chrom_len} --vcftype={config['vcftype']} \
    --non-major-cutoff={config['non_major_cutoff']}; then
    echo "FAILED: ${{peak}}" >&2
    fail=1
  fi
done < pheno.tsv
if [ ${{fail}} -ne 0 ]; then echo 'one or more peaks failed, not writing partial output' >&2; exit 1; fi
first=1
while IFS=$'\\t' read -r peak pheno; do
  if [ ${{first}} -eq 1 ]; then
    printf 'peak\\t' > combined.tsv
    head -n 1 results/${{peak}}.tsv >> combined.tsv
    first=0
  fi
  tail -n +2 results/${{peak}}.tsv | awk -v p=${{peak}} 'BEGIN{{OFS="\\t"}}{{print p, $0}}' >> combined.tsv
done < pheno.tsv
cp combined.tsv {assoc_job.association_results['tsv']}
""",
                )
                b.write_output(assoc_job.association_results, f'{results_dir}/{out_name}')
                manage_concurrency_for_job(assoc_job)
                submitted += 1

            if submitted:
                print(f'{cell_type} {chrom}: submitting {submitted} chunks ({len(entries)} peaks)')
                b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
