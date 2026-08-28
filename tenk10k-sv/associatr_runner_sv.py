#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments,too-many-locals
"""
SV-specific associaTR runner.

Difference vs. str/associatr/associatr_runner.py:
- Uses INFO/END from the EH-style SV VCF (produced by
  tenk10k-sv/sv_vcf_for_associatr.py) to include any SV whose interval
  [POS, END] overlaps the gene cis window by >=1 bp, rather than only variants
  whose POS falls inside the window. This is done per gene with a bcftools
  pre-filter before associaTR is invoked:
      bcftools view -r {chrom}:1-{cis_end} -i 'INFO/END>={cis_start}' ...
  associaTR then walks every record left in the pre-filtered VCF.
- Assumes the cis window bed files were produced with a 1 Mb window
  (see tenk10k-sv/get_cis_numpy_files_sv.toml).

ZERO-VARIANCE PRE-FILTER
------------------------
associaTR standardises the summed genotype as (x - mean(x)) / std(x). A locus with zero dosage
variance among the samples actually used in the regression makes std 0, every value NaN, and
OLS(missing='drop') then discards every row, so statsmodels raises

    ValueError: zero-size array to reduction operation maximum which has no identity

at trtools/associaTR/associaTR.py:280 -- killing the whole gene rather than just that locus.

associaTR's own filters cannot catch it because they work on allele LENGTH frequencies:
sv_vcf_for_associatr.py encodes a multiallelic CNV as GT "CN/0", so a locus where every used
sample has the same copy number still presents as two alleles at ~50% each. It clears the
non-major-allele cutoff and reaches the regression with a constant dosage.

This bites SVs specifically, and rare cell types hardest, because sv_vcf_for_associatr.py drops
monomorphic sites using ALL cohort samples -- a variant that varies across the full cohort but is
constant within one cell type's subset survives into the VCF.

So PREFILTER_SCRIPT below runs inside the associaTR job (the trtools image already has cyvcf2 and
numpy, and the filtered VCF is small) and drops exactly those loci. It imports trtools' own sample
selection, genotype loading and locus filters rather than reimplementing them, and keeps any locus
associaTR would have filtered on its own, so the output TSV is unchanged apart from the omitted
crashers.

Because the pre-filter rewrites the VCF uncompressed and unindexed, --region is no longer passed.
That is equivalent here: bcftools has already restricted the file to the cis window, and the only
thing --region does beyond an index seek is skip records with POS < region_start, which with the
':1-' start never fired.

 analysis-runner --dataset tenk10k --config associatr_runner_sv.toml \
    --description "run associatr on SVs" \
    --access-level test \
    --memory "8G" \
    --output-dir "tenk10k-sv/sv/sensitivity_analysis/tob_n935/1pc" \
    python3 associatr_runner_sv.py
"""

import json

import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path


def gene_cis_window_file_reader(file_path):
    cis_window = pd.read_csv(file_path, sep='\t', header=None)
    cis_window.columns = ['chrom', 'start', 'end']
    return cis_window['chrom'][0], int(cis_window['start'][0]), int(cis_window['end'][0])


# Written into the associaTR job and run against the bcftools-filtered VCF before associaTR itself.
# See the "ZERO-VARIANCE PRE-FILTER" note in the module docstring for why this is needed.
PREFILTER_SCRIPT = r'''
"""Drop the loci that would crash associaTR; pass everything else through untouched.

associaTR standardises the summed genotype as (x - mean(x)) / std(x). If a locus has zero dosage
variance among the samples actually used in the regression, std is 0, every value becomes NaN,
OLS(missing='drop') discards every row, and statsmodels raises

    ValueError: zero-size array to reduction operation maximum which has no identity

which kills the whole gene, not just that locus.

associaTR's own locus filters cannot catch this because they operate on allele LENGTH frequencies.
sv_vcf_for_associatr.py encodes a multiallelic CNV as GT "CN/0", so even when every used sample has
the same copy number the locus presents as two alleles at ~50% each and sails through the
non-major-allele cutoff -- while the summed dosage it hands the regression is constant.

Sample selection, genotype loading and the locus filters are all imported from trtools rather than
reimplemented, so this cannot drift from what associaTR itself does. A locus associaTR would have
filtered anyway is KEPT, so the output TSV still carries its filtered row and is unchanged apart
from the omitted crashers.
"""
import sys

import cyvcf2
import numpy as np

from trtools.associaTR.associaTR import _merge_arrays
from trtools.associaTR.load_and_filter_genotypes import clean_len_alleles
from trtools.utils import tr_harmonizer as trh

vcf_in, npy_path, vcf_out = sys.argv[1], sys.argv[2], sys.argv[3]
non_major_cutoff = float(sys.argv[4])

vcf = cyvcf2.VCF(vcf_in)
all_samples = vcf.samples

# exactly associaTR.py's covariate merge and sample filter (no --sample-list is passed)
covars = np.load(npy_path)
covars = _merge_arrays(np.array(all_samples, dtype=float).reshape(-1, 1), covars)
sample_filter = ~np.any(np.isnan(covars), axis=1)
n_covars = covars.shape[1]

vcftype = trh.InferVCFType(vcf, 'eh')
writer = cyvcf2.Writer(vcf_out, vcf)

n_kept = 0
n_dropped = 0
for record in vcf:
    trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype=vcftype)
    curr_samples = sample_filter & trrecord.GetCalledSamples()
    n_samples = int(np.sum(curr_samples))

    # would associaTR filter this locus before it ever reaches the regression?
    allele_frequency = clean_len_alleles(trrecord.GetAlleleFreqs(curr_samples))
    if len(allele_frequency) <= 1:
        already_filtered = True
    else:
        af = list(allele_frequency.values())
        af.pop(int(np.argmax(af)))
        already_filtered = bool(np.sum(af) * n_samples * 2 < non_major_cutoff)
    already_filtered = already_filtered or n_covars >= n_samples

    if not already_filtered:
        summed_gts = np.sum(trrecord.GetLengthGenotypes()[curr_samples, :-1], axis=1)
        if np.std(summed_gts) == 0:
            print(
                'DROP zero-variance locus {}:{} {} (dosage {} in all {} used samples)'.format(
                    record.CHROM, record.POS, record.ID, summed_gts[0], n_samples,
                ),
            )
            n_dropped += 1
            continue

    writer.write_record(record)
    n_kept += 1

writer.close()
print(
    'prefilter: kept {}, dropped {} zero-variance loci ({} VCF samples, {} usable)'.format(
        n_kept, n_dropped, len(all_samples), int(np.sum(sample_filter)),
    ),
)
'''


def main():
    b = get_batch(name='Run associatr on SVs')

    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        if len(_dependent_jobs) >= get_config()['associatr']['max_parallel_jobs']:
            job.depends_on(_dependent_jobs[-get_config()['associatr']['max_parallel_jobs']])
        _dependent_jobs.append(job)

    for celltype in get_config()['associatr']['celltypes'].split(','):
        for chromosome in get_config()['associatr']['chromosomes'].split(','):
            input_dir = get_config()['associatr']['vcf_file_dir']
            vcf_file_path = f'{input_dir}/hail_filtered_{chromosome}.vcf.bgz'
            variant_vcf = b.read_input_group(
                base=vcf_file_path,
                tbi=vcf_file_path + '.tbi',
            )

            gene_list_dir = get_config()['associatr']['gene_list_dir']
            with to_path(f'{gene_list_dir}/{celltype}/{chromosome}_{celltype}_gene_list.json').open('r') as file:
                pseudobulk_gene_names = json.load(file)

            for gene in pseudobulk_gene_names:
                cis_window_dir = get_config()['associatr']['cis_window_dir']
                cis_window_size = get_config()['associatr']['cis_window_size']
                version = get_config()['associatr']['version']

                if to_path(
                    output_path(
                        f'results/{version}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.tsv',
                        'analysis',
                    ),
                ).exists():
                    continue

                gene_cis_window_file = f'{cis_window_dir}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp.bed'
                cis_chrom, cis_start, cis_end = gene_cis_window_file_reader(gene_cis_window_file)

                pheno_cov_numpy_dir = get_config()['associatr']['pheno_cov_numpy_dir']
                gene_pheno_cov = b.read_input(
                    f'{pheno_cov_numpy_dir}/{celltype}/{chromosome}/{gene}_pheno_cov.npy',
                )

                # Job 1: bcftools pre-filter to the SVs whose [POS, INFO/END]
                # intersects [cis_start, cis_end] by >=1 bp. Runs on the
                # bcftools image because the trtools image does not ship
                # bcftools/htslib. Output is a bgzip+tabix pair carried through
                # the batch as a resource group to the associaTR job below.
                filter_job = b.new_job(name=f'Filter SV VCF for {gene} [{celltype};{chromosome}]')
                if get_config()['associatr']['always_run']:
                    filter_job.always_run()
                filter_job.image(get_config()['images']['bcftools'])
                filter_job.storage(get_config()['associatr']['job_storage'])
                filter_job.cpu(get_config()['associatr']['job_cpu'])
                filter_job.declare_resource_group(
                    filtered_vcf={
                        'vcf.bgz': '{root}.vcf.bgz',
                        'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
                    },
                )
                filter_job.command(
                    f"bcftools view -r {cis_chrom}:1-{cis_end} "
                    f"-i 'INFO/END>={cis_start}' "
                    f"{variant_vcf.base} -Oz -o {filter_job.filtered_vcf['vcf.bgz']} && "
                    f"tabix -f -p vcf {filter_job.filtered_vcf['vcf.bgz']}",
                )

                # Job 2: associaTR on the per-gene filtered VCF, preceded by the
                # zero-variance pre-filter (see PREFILTER_SCRIPT). The pre-filter
                # rewrites the VCF uncompressed and unindexed, so --region is
                # dropped: it is equivalent here anyway, because bcftools has
                # already restricted the file to the cis window and the only thing
                # --region does beyond an index seek is skip records starting
                # before region_start -- which with ':1-' never fired.
                non_major_cutoff = get_config()['associatr'].get('non_major_cutoff', 20)
                associatr_job = b.new_job(name=f'Run associatr on {gene} [{celltype};{chromosome}] SV')
                if get_config()['associatr']['always_run']:
                    associatr_job.always_run()
                associatr_job.image(get_config()['images']['trtools'])
                associatr_job.storage(get_config()['associatr']['job_storage'])
                associatr_job.cpu(get_config()['associatr']['job_cpu'])
                associatr_job.declare_resource_group(
                    association_results={'tsv': '{root}.tsv'},
                )
                associatr_job.command(
                    f'cat > $BATCH_TMPDIR/prefilter.py <<\'PREFILTER_EOF\'\n'
                    f'{PREFILTER_SCRIPT}\n'
                    f'PREFILTER_EOF',
                )
                associatr_job.command(
                    f'python3 $BATCH_TMPDIR/prefilter.py '
                    f"{filter_job.filtered_vcf['vcf.bgz']} {gene_pheno_cov} "
                    f'$BATCH_TMPDIR/prefiltered.vcf {non_major_cutoff}',
                )
                associatr_job.command(
                    f"associaTR {associatr_job.association_results['tsv']} "
                    f'$BATCH_TMPDIR/prefiltered.vcf '
                    f"{celltype}_{chromosome}_{gene} {gene_pheno_cov} "
                    f'--non-major-cutoff {non_major_cutoff} --vcftype=eh',
                )
                b.write_output(
                    associatr_job.association_results,
                    output_path(
                        f'results/{version}/{celltype}/{chromosome}/{gene}_{cis_window_size}bp',
                        'analysis',
                    ),
                )
                manage_concurrency_for_job(associatr_job)
    b.run(wait=False)


if __name__ == '__main__':
    main()
