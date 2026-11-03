#!/usr/bin/env python3

"""
Drop F <= -0.5 SVs from the SV eQTL meta-analysis results.

A federated pair is removed if EITHER side failed the inbreeding-coefficient filter in its own
cohort: the TOB VARID (`tob_vid`) is listed in the TOB excluded table, or the BioHEART VARID
(`bioheart_vid`) is listed in the BioHEART one. The excluded tables are the
*.F_le_minus0.5.excluded.tsv files written alongside the F-filtered VCFs (`ID` column = VCF VARID).

Input layout (output of meta_runner_sv.py):
    {meta-results-dir}/{cell_type}/{chrom}/{gene}_{cis_window}bp_meta_results.tsv

Output, same layout under the analysis output dir:
    {output-dir}/meta_results/{cell_type}/{chrom}/{gene}_{cis_window}bp_meta_results.tsv

Files are filtered line by line, so surviving rows are byte-identical to the input (and keep the
pval_meta_fixed sort order). A gene whose rows are all excluded gets no output file, matching
meta_runner_sv.py, which writes nothing for genes with no poolable pair. gene_level_pvals and
fdr_qvals must be regenerated from the new meta_results.

analysis-runner --dataset "tenk10k" --description "SV meta-analysis F filter" --access-level "test" \
    --output-dir "tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v1_F_filtered" \
    --image 'australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:189f33cce2381aa8623acd2f39006455deae2a6f-hail-0928b141d1dd3ca4de530d50943d2673305ea7bb' \
    python3 meta_results_F_filter_sv.py \
    --meta-results-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/meta_analysis/bioheart_n968_and_tob_n935/7pc/v1/meta_results \
    --tob-excluded=gs://cpg-tenk10k-test/sv/input_files/F_filtered_vcfs/tob_common_maf_gte_1pct.F_le_minus0.5.excluded.tsv \
    --bioheart-excluded=gs://cpg-tenk10k-test/sv/input_files/F_filtered_vcfs/bioheart_common_maf_gte_1pct.F_le_minus0.5.excluded.tsv \
    --cell-types=ASDC,B_intermediate,B_memory,B_naive,CD14_Mono,CD16_Mono,CD4_CTL,CD4_Naive,CD4_Proliferating,CD4_TCM,CD4_TEM,CD8_Naive,CD8_Proliferating,CD8_TCM,CD8_TEM,HSPC,ILC,MAIT,NK,NK_CD56bright,NK_Proliferating,Plasmablast,Treg,cDC1,cDC2,dnT,gdT,pDC \
    --chromosomes=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
"""

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path


def _load_excluded(path):
    """excluded.tsv -> set of upper-cased VARIDs (case-insensitive match, as in meta_runner_sv.py)."""
    with to_path(path).open() as handle:
        ids = pd.read_csv(handle, sep='\t', usecols=['ID'])['ID']
    return set(ids.astype(str).str.strip().str.upper())


def filter_chrom(meta_results_dir, tob_excluded, bioheart_excluded, cell_type, chromosome):
    """Filter every gene meta-results file for one cell type / chromosome."""
    from cpg_utils.hail_batch import output_path

    tob_ex = _load_excluded(tob_excluded)
    bh_ex = _load_excluded(bioheart_excluded)

    in_files = sorted(to_path(f'{meta_results_dir}/{cell_type}/{chromosome}').glob('*_meta_results.tsv'))
    print(f'{cell_type} {chromosome}: {len(in_files):,} gene files')

    n_rows_in = n_rows_out = n_files_out = n_emptied = n_skipped = 0
    for in_file in in_files:
        out_path = output_path(f'meta_results/{cell_type}/{chromosome}/{in_file.name}', 'analysis')
        if to_path(out_path).exists():
            n_skipped += 1
            continue

        lines = in_file.read_text().splitlines(keepends=True)
        header = lines[0].rstrip('\n').split('\t')
        i_tob, i_bh = header.index('tob_vid'), header.index('bioheart_vid')

        kept = []
        for line in lines[1:]:
            fields = line.rstrip('\n').split('\t')
            if fields[i_tob].upper() in tob_ex or fields[i_bh].upper() in bh_ex:
                continue
            kept.append(line)
        n_rows_in += len(lines) - 1
        n_rows_out += len(kept)

        if not kept:
            n_emptied += 1
            continue
        with to_path(out_path).open('w') as handle:
            handle.write(lines[0])
            handle.writelines(kept)
        n_files_out += 1

    print(
        f'  rows {n_rows_in:,} -> {n_rows_out:,} ({n_rows_in - n_rows_out:,} dropped) | '
        f'wrote {n_files_out:,} files, {n_emptied:,} genes emptied, {n_skipped:,} already present',
    )


@click.command()
@click.option('--meta-results-dir', required=True, help='meta_results dir written by meta_runner_sv.py')
@click.option('--tob-excluded', required=True, help='TOB *.F_le_minus0.5.excluded.tsv (ID column)')
@click.option('--bioheart-excluded', required=True, help='BioHEART *.F_le_minus0.5.excluded.tsv (ID column)')
@click.option('--cell-types', required=True, help='comma-separated cell types')
@click.option('--chromosomes', required=True, help='comma-separated chromosomes')
@click.option('--max-parallel-jobs', default=500, type=int, help='concurrency cap')
@click.option('--always-run', is_flag=True, help='set jobs to always run')
def main(meta_results_dir, tob_excluded, bioheart_excluded, cell_types, chromosomes, max_parallel_jobs, always_run):
    """Fan out one job per cell type / chromosome."""
    from cpg_utils.hail_batch import get_batch

    # fail fast on bad paths / missing ID column before spending any compute
    print(f'TOB excluded: {len(_load_excluded(tob_excluded)):,} | BioHEART excluded: {len(_load_excluded(bioheart_excluded)):,}')

    b = get_batch(name='SV meta-analysis F filter')

    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """Limit how many jobs run at once."""
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    for cell_type in cell_types.split(','):
        for chromosome in chromosomes.split(','):
            j = b.new_python_job(name=f'F_filter_{cell_type}_{chromosome}')
            if always_run:
                j.always_run()
            j.call(filter_chrom, meta_results_dir, tob_excluded, bioheart_excluded, cell_type, chromosome)
            manage_concurrency_for_job(j)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
