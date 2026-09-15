#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
Export AGGREGATED genotype-vs-expression tables for exemplar loci (deliverable 8).

This is the shareable replacement for plotters/sv_genotype_vs_residual_expression.ipynb. That
notebook plots one point per donor, which is individual-level data. This script emits one row per
DOSAGE BIN:

    cohort, cell_type, gene_name, varid, dosage, n_donors, mean_residual_expr, se_residual_expr

so the same effect plot can be drawn -- means with error bars instead of a point cloud -- without
any individual's genotype or expression leaving the environment. Bins holding fewer than
--min-donors donors are dropped entirely (not merged), and the count of suppressed bins is reported
so the suppression itself is auditable.

WHAT "RESIDUAL EXPRESSION" MEANS HERE
-------------------------------------
The phenotype is residualised on the covariates only, and binned against RAW dosage, because a raw
x-axis (0/1/2 alt alleles, or copy number) is what makes the plot readable.

That is NOT identical to associaTR's model, which fits phenotype on dosage and covariates jointly.
By Frisch-Waugh-Lovell the joint slope is recovered only when BOTH phenotype and dosage are
residualised on the covariates. So the script also computes that doubly-residualised slope per locus
and writes it as `fwl_slope`, which should reproduce the associaTR `coeff` for the same cohort to
numerical precision. If it does not, the sample selection has drifted and the binned means should
not be trusted -- check that before using the output.

SAMPLE SELECTION
----------------
Imported from trtools rather than reimplemented, exactly as the zero-variance pre-filter in
associatr_runner_sv.py does, so which donors enter a bin cannot drift from which donors entered the
regression.

DOSAGE ENCODING
---------------
Inherited from sv_vcf_for_associatr.py: biallelic records give alt allele count (0/1/2), while
FILTER=MULTIALLELIC CNV records give copy number (GT "CN/0" summing to CN). So `dosage` is per-allele
for the former and per-copy for the latter, and a CNV's reference bin sits at 2, not 0. Do not put
the two on a shared axis.

analysis-runner --dataset "tenk10k" --description "export dosage bins" --access-level "test" \
    --output-dir "tenk10k-sv/sv/export/dosage_bins/v1" \
    --image 'australia-southeast1-docker.pkg.dev/cpg-common/images/trtools:6.0.1' \
    python3 extract_locus_dosage_bins.py \
    --requests=gs://cpg-tenk10k-test/sv/export/dosage_bin_requests.tsv \
    --vcf-dir=gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/bioheart_common_chr_specific \
    --pheno-cov-dir=gs://cpg-tenk10k-test/sv/input_files/bioheart_n968/7pc/pheno_cov_numpy/v1 \
    --cohort=bioheart

REQUESTS FILE
-------------
Tab- or comma-separated, one row per locus to export. Columns:

    cell_type   chr     gene_name         varid                           pos
    NK          chr1    ENSG00000187010   all_batches.chr1..._DEL_1_725   25272393

`varid` and `pos` are the per-cohort values (bioheart_vid/bioheart_pos or tob_vid/tob_pos from the
meta-analysis output) -- NOT federate_id, which has no entry in either cohort's VCF. `pos` is used
for an indexed region seek, so it must match the VCF record.
"""

import json

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

# Written into the job and run against the cohort's EH-style SV VCF. Kept as a shipped script (rather
# than a PythonJob) because it needs cyvcf2 + trtools, which live in the trtools image, not the
# driver image -- the same reason associatr_runner_sv.py ships its pre-filter this way.
DOSAGE_BIN_SCRIPT = r'''
"""Aggregate covariate-residualised expression by raw SV dosage. Writes no per-donor value."""
import json
import sys

import cyvcf2
import numpy as np

from trtools.associaTR.associaTR import _merge_arrays
from trtools.utils import tr_harmonizer as trh

vcf_path, npy_path, requests_json, out_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
min_donors = int(sys.argv[5])
cohort, cell_type, gene = sys.argv[6], sys.argv[7], sys.argv[8]

requests = json.loads(open(requests_json).read())  # [{"varid": ..., "pos": ...}, ...]

vcf = cyvcf2.VCF(vcf_path)
all_samples = vcf.samples

# associaTR's covariate merge and sample filter, verbatim (no --sample-list is passed).
# Column layout of the .npy, set by str/associatr/get_cis_numpy_files.py:
#   0 = sample id (numeric), 1 = INT-transformed expression, 2+ = covariates
covars = np.load(npy_path)
covars = _merge_arrays(np.array(all_samples, dtype=float).reshape(-1, 1), covars)
sample_filter = ~np.any(np.isnan(covars), axis=1)

vcftype = trh.InferVCFType(vcf, 'eh')

wanted = {str(r['varid']).upper(): int(r['pos']) for r in requests}
rows = []
n_suppressed = 0

for varid, pos in wanted.items():
    # indexed seek; a 2 bp window absorbs any off-by-one between the results table and the VCF
    found = None
    for record in vcf(f'{record_chrom_placeholder}:{max(1, pos - 1)}-{pos + 1}'):
        if str(record.ID).upper() == varid:
            found = record
            break
    if found is None:
        print(f'  !! {varid} not found at pos {pos}; skipping')
        continue

    trrecord = trh.HarmonizeRecord(vcfrecord=found, vcftype=vcftype)
    curr = sample_filter & trrecord.GetCalledSamples()
    n_used = int(np.sum(curr))
    if n_used == 0:
        print(f'  !! {varid} has no usable samples; skipping')
        continue

    dosage = np.sum(trrecord.GetLengthGenotypes()[curr, :-1], axis=1)
    pheno = covars[curr, 1]
    cov = covars[curr, 2:]

    if np.std(dosage) == 0:
        print(f'  !! {varid} has zero dosage variance in {cell_type}; skipping')
        continue

    # residualise the phenotype on [intercept, covariates]
    design = np.column_stack([np.ones(n_used), cov])
    beta, *_ = np.linalg.lstsq(design, pheno, rcond=None)
    resid_pheno = pheno - design @ beta

    # Frisch-Waugh-Lovell check: residualise dosage on the SAME design, then regress. This slope
    # must equal associaTR's coeff for this cohort/locus/gene. The binned means below use RAW
    # dosage, so their apparent slope can differ slightly -- this number is the ground truth.
    beta_d, *_ = np.linalg.lstsq(design, dosage, rcond=None)
    resid_dosage = dosage - design @ beta_d
    denom = float(np.sum(resid_dosage ** 2))
    fwl_slope = float(np.sum(resid_dosage * resid_pheno) / denom) if denom > 0 else float('nan')

    for d in np.unique(dosage):
        mask = dosage == d
        n = int(np.sum(mask))
        if n < min_donors:
            n_suppressed += 1
            continue
        vals = resid_pheno[mask]
        # SE of the mean; undefined for a single donor, but n >= min_donors >= 2 by validation
        se = float(np.std(vals, ddof=1) / np.sqrt(n))
        rows.append({
            'cohort': cohort,
            'cell_type': cell_type,
            'gene_name': gene,
            'varid': varid,
            'dosage': float(d),
            'n_donors': n,
            'mean_residual_expr': float(np.mean(vals)),
            'se_residual_expr': se,
            'n_donors_total': n_used,
            'fwl_slope': fwl_slope,
        })

with open(out_path, 'w') as f:
    cols = ['cohort', 'cell_type', 'gene_name', 'varid', 'dosage', 'n_donors',
            'mean_residual_expr', 'se_residual_expr', 'n_donors_total', 'fwl_slope']
    f.write('\t'.join(cols) + '\n')
    for r in rows:
        f.write('\t'.join(str(r[c]) for c in cols) + '\n')

print(f'wrote {len(rows)} dosage bins for {cell_type}/{gene} '
      f'({n_suppressed} bins suppressed with n < {min_donors})')
'''


def read_requests(path):
    """Requests file -> DataFrame, validated."""
    with to_path(path).open() as handle:
        sep = '\t' if str(path).endswith(('.tsv', '.txt')) else ','
        req = pd.read_csv(handle, sep=sep)

    required = {'cell_type', 'chr', 'gene_name', 'varid', 'pos'}
    missing = required - set(req.columns)
    if missing:
        raise ValueError(f'{path} is missing column(s): {sorted(missing)}')

    req['chr'] = req['chr'].apply(lambda c: c if str(c).startswith('chr') else f'chr{c}')
    print(f'{len(req):,} locus requests across {req.cell_type.nunique()} cell types, '
          f'{req.gene_name.nunique()} genes, {req.varid.nunique()} SVs')
    return req


@click.command()
@click.option('--requests', 'requests_path', required=True, help='GCS path to the requests file')
@click.option('--vcf-dir', required=True, help='GCS dir of per-chromosome EH-style SV VCFs (.vcf.bgz + .tbi)')
@click.option('--pheno-cov-dir', required=True, help='GCS dir of {cell_type}/{chrom}/{gene}_pheno_cov.npy')
@click.option('--cohort', required=True, help='Cohort label written into the output')
@click.option('--min-donors', default=5, type=int, help='drop any dosage bin with fewer donors')
@click.option('--job-storage', default='5G')
@click.option('--job-cpu', default=1.0, type=float)
def main(requests_path, vcf_dir, pheno_cov_dir, cohort, min_donors, job_storage, job_cpu):
    """Export aggregated dosage-binned expression for a hand-picked set of loci."""
    if min_donors < 2:
        raise click.UsageError('--min-donors must be >= 2; the SE of a single-donor bin is undefined')

    req = read_requests(requests_path)
    b = get_batch(name=f'Dosage bin export ({cohort})')
    outputs = []

    # read each chromosome's VCF once and share the resource across every job that needs it,
    # rather than re-declaring it per gene
    vcfs = {
        chrom: b.read_input_group(
            **{'vcf.bgz': f'{vcf_dir}/{chrom}.vcf.bgz', 'vcf.bgz.tbi': f'{vcf_dir}/{chrom}.vcf.bgz.tbi'},
        )
        for chrom in sorted(req['chr'].unique())
    }

    # one job per (cell_type, chrom, gene): each needs its own pheno_cov .npy, and that file is the
    # unit the association itself was run on
    for (cell_type, chrom, gene), grp in req.groupby(['cell_type', 'chr', 'gene_name']):
        vcf = vcfs[chrom]
        npy = b.read_input(f'{pheno_cov_dir}/{cell_type}/{chrom}/{gene}_pheno_cov.npy')

        j = b.new_job(name=f'dosage bins {gene} [{cell_type};{chrom}]')
        j.image(get_config()['images']['trtools'])
        j.storage(job_storage)
        j.cpu(job_cpu)

        payload = json.dumps([{'varid': v, 'pos': int(p)} for v, p in zip(grp.varid, grp.pos)])
        # chromosome is fixed per job, so the region string is templated in here rather than
        # threaded through argv
        script = DOSAGE_BIN_SCRIPT.replace('{record_chrom_placeholder}', chrom)

        j.command(f"cat > $BATCH_TMPDIR/bins.py <<'BINS_EOF'\n{script}\nBINS_EOF")
        j.command(f"cat > $BATCH_TMPDIR/req.json <<'REQ_EOF'\n{payload}\nREQ_EOF")
        j.command(
            f'python3 $BATCH_TMPDIR/bins.py '
            f"{vcf['vcf.bgz']} {npy} $BATCH_TMPDIR/req.json {j.out} "
            f'{min_donors} {cohort} {cell_type} {gene}',
        )
        outputs.append((f'{cohort}/{cell_type}/{chrom}/{gene}_dosage_bins.tsv', j.out))

    for path, resource in outputs:
        b.write_output(resource, output_path(path, 'analysis'))

    print(f'submitting {len(outputs)} jobs; bins with < {min_donors} donors are dropped')
    print('Output rows are group means and counts only -- no per-donor values.')
    b.run(wait=False)


if __name__ == '__main__':
    main()
