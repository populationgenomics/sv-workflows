#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
Export SITE-LEVEL tables from a GATK-SV VCF (deliverables 1 and 2).

Produces two files, neither of which contains per-sample data:

    sv_sites.tsv.gz      one row per SV: coordinates, SVLEN, SVTYPE, FILTER, allele counts,
                         calling evidence, and the svtk PREDICTED_* consequence annotations
                         when the VCF carries them
    sv_site_qc.tsv.gz    one row per SV: recomputed AC/AN/AF plus F_MISSING and HWE

Everything here is an aggregate over samples (counts, frequencies, test statistics). No genotype,
dosage or copy-number value for any individual is written. That is the point of the script: it is
the export that can leave a controlled-access environment.

WHY THIS FILE IS NEEDED
-----------------------
sv_vcf_for_associatr.py keeps only SVTYPE and END out of the GATK-SV INFO block, packing them into
the mock ExpansionHunter `RU` field. Every other annotation -- SVLEN above all -- is dropped before
association, and never appears in the associaTR output or the meta-analysis results. So downstream
there is no way to ask:

  - how long is this SV?             END - POS is ~0 for insertions; only SVLEN answers it
  - does it overlap the gene?        needs the interval, or the PREDICTED_* annotations
  - was it called from depth alone?  needs ALGORITHMS / EVIDENCE

sv_sites.tsv is the join-back table: key on ID (the VARID that associaTR carries through the `motif`
column, upper-cased) and every one of those questions becomes answerable against results that are
already computed.

PREDICTED_* FIELDS
------------------
These come from GATK-SV's AnnotateVcf / `svtk annotate` step. If that step was not run the VCF has
none of them, and the script simply omits those columns rather than failing -- the header is
inspected at runtime and the bcftools format string is built from what is actually present. Gene
overlap can then only be derived from the interval, which is coarser (no exon/UTR/promoter
resolution), so running AnnotateVcf first is worth the trouble if it is available.

F_MISSING AND HWE ON MULTIALLELIC CNVs
--------------------------------------
`bcftools +fill-tags` derives its tags from GT. GATK-SV does not emit meaningful GT for records
FILTERed as MULTIALLELIC -- copy number lives in FORMAT/CN instead -- so AC/AN/AF/F_MISSING/HWE are
not interpretable for those rows. FILTER is exported alongside so they can be excluded downstream;
use the per-site dosage spectrum (the `allele_frequency` JSON already present in the associaTR
output) for multiallelic CNVs instead.

analysis-runner --dataset "tenk10k-sv" --description "export SV site tables" --access-level "test" \
    --output-dir "pub-analysis/export/sv_sites/v1" \
    python3 extract_sv_sites.py \
    --vcf-path=gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/bioheart_common_maf_gte_1pct.vcf.gz \
    --cohort=bioheart
"""

import click

from cpg_utils.config import get_config, output_path
from cpg_utils.hail_batch import get_batch

# INFO fields pulled unconditionally. Everything here is either mandatory in a VCF or emitted by
# GATK-SV on every record; a missing value comes back as '.' rather than failing the query.
CORE_INFO_FIELDS = [
    ('ID', '%ID'),
    ('chrom', '%CHROM'),
    ('pos', '%POS'),
    ('end', '%INFO/END'),
    ('svlen', '%INFO/SVLEN'),
    ('svtype', '%INFO/SVTYPE'),
    ('filter', '%FILTER'),
    ('ac', '%INFO/AC'),
    ('an', '%INFO/AN'),
    ('af', '%INFO/AF'),
    ('algorithms', '%INFO/ALGORITHMS'),
    ('evidence', '%INFO/EVIDENCE'),
]

# svtk/AnnotateVcf consequence annotations. Probed against the VCF header at runtime; any that the
# VCF does not declare are left out of the query rather than emitted as '.' columns, so the output
# says honestly which annotations exist.
OPTIONAL_PREDICTED_FIELDS = [
    'PREDICTED_LOF',
    'PREDICTED_COPY_GAIN',
    'PREDICTED_DUP_PARTIAL',
    'PREDICTED_PARTIAL_EXON_DUP',
    'PREDICTED_INTRAGENIC_EXON_DUP',
    'PREDICTED_TSS_DUP',
    'PREDICTED_INTRONIC',
    'PREDICTED_UTR',
    'PREDICTED_PROMOTER',
    'PREDICTED_BREAKEND_EXONIC',
    'PREDICTED_INV_SPAN',
    'PREDICTED_MSV_EXON_OVERLAP',
    'PREDICTED_NEAREST_TSS',
    'PREDICTED_INTERGENIC',
    'PREDICTED_NONCODING_SPAN',
    'PREDICTED_NONCODING_BREAKPOINT',
]

QC_FIELDS = [
    ('ID', '%ID'),
    ('filter', '%FILTER'),
    ('ac_recomputed', '%INFO/AC'),
    ('an_recomputed', '%INFO/AN'),
    ('af_recomputed', '%INFO/AF'),
    ('f_missing', '%INFO/F_MISSING'),
    ('hwe', '%INFO/HWE'),
]


def _query_expr(fields):
    """[(name, spec), ...] -> (tab-joined header line, bcftools -f body WITHOUT a trailing newline).

    The newline is left off so callers can append further fields (the PREDICTED_* block) before
    terminating the record.
    """
    header = '\t'.join(name for name, _ in fields)
    body = '\\t'.join(spec for _, spec in fields)
    return header, body


def build_command(cohort: str, sites_out: str, qc_out: str) -> str:
    """Bash run inside the bcftools job.

    The PREDICTED_* probe happens here rather than in the driver because the driver would have to
    stream the VCF header out of GCS to do it, and the job already has the file localised.
    """
    core_header, core_body = _query_expr(CORE_INFO_FIELDS)
    qc_header, qc_body = _query_expr(QC_FIELDS)
    predicted_list = ' '.join(OPTIONAL_PREDICTED_FIELDS)

    return f"""
set -euo pipefail

echo "=== probing header for PREDICTED_* annotations ==="
PRESENT=""
HEADER_EXTRA=""
for f in {predicted_list}; do
  if bcftools view -h "$VCF" | grep -q "ID=${{f}},"; then
    PRESENT="${{PRESENT}}\\\\t%INFO/${{f}}"
    HEADER_EXTRA="${{HEADER_EXTRA}}\\t${{f}}"
    echo "  found $f"
  fi
done
if [ -z "$PRESENT" ]; then
  echo "  none found -- VCF has not been through AnnotateVcf/svtk annotate."
  echo "  Gene overlap will have to be derived from POS/END/SVLEN downstream."
fi

echo "=== 1/2 sv_sites ==="
# %b so the \\t separators accumulated in HEADER_EXTRA become real tabs
printf '{core_header}%b\\n' "$HEADER_EXTRA" > sites.tsv
bcftools query -f '{core_body}'"$PRESENT"'\\n' "$VCF" >> sites.tsv
gzip -c sites.tsv > {sites_out}

echo "=== 2/2 sv_site_qc ==="
# +fill-tags RECOMPUTES AC/AN/AF from GT, so these columns are deliberately named *_recomputed and
# will disagree with the INFO/AC carried in sv_sites for any record whose GT is not meaningful
# (notably FILTER=MULTIALLELIC). That disagreement is informative, so both are kept.
printf '{qc_header}\\n' > qc.tsv
bcftools +fill-tags "$VCF" -Ou -- -t AC,AN,AF,F_MISSING,HWE \\
  | bcftools query -f '{qc_body}\\n' >> qc.tsv
gzip -c qc.tsv > {qc_out}

echo "=== done: {cohort} ==="
wc -l sites.tsv qc.tsv
echo "Both outputs are site-level aggregates; no per-sample values are written."
"""


@click.command()
@click.option('--vcf-path', required=True, help='GCS path to the GATK-SV VCF (bgzipped)')
@click.option('--cohort', required=True, help='Cohort label, used in the output filenames')
@click.option('--job-storage', default='30G')
@click.option('--job-cpu', default=2)
def main(vcf_path, cohort, job_storage, job_cpu):
    """Export site-level SV tables for sharing outside the analysis environment."""
    b = get_batch(name=f'SV site export ({cohort})')

    j = b.new_job(name=f'bcftools site query [{cohort}]')
    j.image(get_config()['images']['bcftools'])
    j.storage(job_storage)
    j.cpu(job_cpu)
    j.declare_resource_group(
        out={
            'sites.tsv.gz': '{root}.sites.tsv.gz',
            'site_qc.tsv.gz': '{root}.site_qc.tsv.gz',
        },
    )

    vcf = b.read_input(vcf_path)
    j.command(f'export VCF={vcf}')
    j.command(build_command(cohort, j.out['sites.tsv.gz'], j.out['site_qc.tsv.gz']))

    b.write_output(j.out, output_path(f'{cohort}', 'analysis'))
    b.run(wait=False)


if __name__ == '__main__':
    main()
