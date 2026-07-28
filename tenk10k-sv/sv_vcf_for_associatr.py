#!/usr/bin/env python3
"""
This script converts an existing GATK-SV VCF file into a
mock ExpansionHunter-style VCF for direct use with associaTR.

Please note ad-hoc changes:
- RU field stores the REF and ALT alleles in the format of 'REF-ALT' (REF/ALT
  here are the SV's symbolic alleles, e.g. 'N-<DEL>', not sequence bases)
- RL field is set to 0 for all loci so that the REF allele is coded as 0
  (https://github.com/gymrek-lab/TRTools/blob/master/trtools/utils/tr_harmonizer.py#L515)
- REF field (INFO) is set to 3 (arbitrary filler value, matches SNP version)
- SVTYPE and END are taken from the original GATK-SV INFO field (not faked),
  since GATK-SV already provides these and they're meaningful for SVs
- GT is extracted from whichever FORMAT slot it occupies (GATK-SV FORMAT is
  typically GT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ:SL), not assumed to be field 0 only
  based on position -- it's located by name in the FORMAT header each time
- Multi-allelic CNVs (FILTER contains MULTIALLELIC) do not carry meaningful
  GT calls in GATK-SV; for these records we read FORMAT/CN (copy number) per
  sample and emit it as "CN/0" in the mock GT slot. associaTR sums the two
  alleles per sample, so effective per-sample dosage = CN. Missing CN is
  written as "./.". The monomorphic-site drop uses CN directly for these
  records. The ALT column for these records is expanded to
  <STR1>,<STR2>,...,<STR{max_cn}> so allele indices in the mock GT resolve
  to a valid allele (otherwise associaTR IndexErrors in GetLengthGenotypes).

analysis-runner --dataset tenk10k-sv --access-level test --output-dir pub-analysis/final-vcf/filtered --description "sv vcf for associatr" \
    sv_vcf_for_associatr.py --vcf-path=gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/tob_common_maf_gte_1pct.vcf.gz \
    --job-storage=10G --job-cpu=8
"""

import gzip
import os
import sys

import click
from cpg_utils import to_path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch


# IDs that the EH-style header block below redefines. Any matching header line
# from the source GATK-SV VCF is skipped on passthrough so the output does not
# contain two definitions for the same INFO/FORMAT id
_REDEFINED_HEADER_IDS = {
    ('INFO', 'END'),
    ('INFO', 'REF'),
    ('INFO', 'REPID'),
    ('INFO', 'RL'),
    ('INFO', 'RU'),
    ('INFO', 'SVTYPE'),
    ('INFO', 'VARID'),
    ('FORMAT', 'GT'),
    ('FORMAT', 'ADFL'),
    ('FORMAT', 'ADIR'),
    ('FORMAT', 'ADSP'),
    ('FORMAT', 'LC'),
    ('FORMAT', 'REPCI'),
    ('FORMAT', 'REPCN'),
    ('FORMAT', 'SO'),
    ('FORMAT', 'QUAL'),
}


def _is_redefined_header(header_line):
    """Return True if header_line declares an INFO/FORMAT id the EH block redefines."""
    for kind in ('INFO', 'FORMAT'):
        prefix = f'##{kind}=<ID='
        if header_line.startswith(prefix):
            hid = header_line[len(prefix):].split(',', 1)[0]
            return (kind, hid) in _REDEFINED_HEADER_IDS
    return False


def _log_drop(reason, chrom, pos, record_id, svtype, extra=''):
    """Emit a single-line drop record to stderr.

    Format: DROP<TAB>reason<TAB>chrom<TAB>pos<TAB>id<TAB>svtype<TAB>extra
    Grep-friendly, one line per dropped variant; also feeds the end-of-run summary.
    """
    print(
        f'DROP\t{reason}\t{chrom}\t{pos}\t{record_id}\t{svtype}\t{extra}',
        file=sys.stderr,
    )


NEW_FORMAT_FIELDS = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n',
    '##FORMAT=<ID=ADFL,Number=1,Type=String,Description="">\n',
    '##FORMAT=<ID=ADIR,Number=1,Type=String,Description="">\n',
    '##FORMAT=<ID=ADSP,Number=1,Type=String,Description="">\n',
    '##FORMAT=<ID=LC,Number=1,Type=Float,Description="">\n',
    '##FORMAT=<ID=REPCI,Number=1,Type=String,Description="">\n',
    '##FORMAT=<ID=REPCN,Number=1,Type=String,Description="">\n',
    '##FORMAT=<ID=SO,Number=1,Type=String,Description="">\n',
    '##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="">\n',
]

NEW_INFO_FIELDS = [
    '##INFO=<ID=END,Number=1,Type=Integer,Description="">\n',
    '##INFO=<ID=REF,Number=1,Type=Integer,Description="">\n',
    '##INFO=<ID=REPID,Number=1,Type=String,Description="">\n',
    '##INFO=<ID=RL,Number=1,Type=Integer,Description="">\n',
    (
        '##INFO=<ID=RU,Number=1,Type=String,Description="Storing the REF/ALT info instead of '
        'Repeat Unit to retain the REF/ALT info in the association output files, which is '
        'crucial where there are multiple variants with the same CHR:POS coordinates '
        '(eg multi allelic loci)">\n'
    ),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">\n',
    '##INFO=<ID=VARID,Number=1,Type=String,Description="">\n',
]


def parse_info_field(info_field):
    """Parse a VCF INFO string into a dict, handling flag-only entries."""
    info = {}
    for entry in info_field.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info[key] = value
        else:
            info[entry] = True
    return info


def reformat_vcf(vcf_file_path, output_file_path):
    vcf_file = to_path(vcf_file_path)
    output_file = to_path(output_file_path)

    drop_counts = {
        'multiallelic_no_cn': 0,
        'monomorphic_no_callable': 0,
        'monomorphic_uniform_dosage': 0,
    }
    n_seen = 0
    n_written = 0

    # Write the temp file gzip-compressed so the uploaded blob's contents match
    # its .vcf.gz extension (inherited from the input basename in main()).
    # A downstream bgzip+tabix step is still required before this VCF can be
    # consumed by associatr_runner_sv.py, which expects a .vcf.bgz + .tbi pair.
    with gzip.open(vcf_file, 'rt') as fin, gzip.open('temporary_gt_file.vcf.gz', 'wt') as fout:
        for line in fin:
            if line.startswith('#CHROM'):
                # Insert the mock ExpansionHunter-style header fields right before
                # the #CHROM line (GATK-SV VCFs don't have a ##hailversion line to
                # key off of, unlike the hail-derived SNP VCFs).
                fout.writelines(NEW_FORMAT_FIELDS)
                fout.writelines(NEW_INFO_FIELDS)
                # Declare a range of <STRn> symbolic alleles up front. MULTIALLELIC
                # CN records emit ALT = <STR1>,<STR2>,...,<STRmax_cn> so allele
                # indices in the mock GT (e.g. "3/0" for CN=3) resolve to a valid
                # allele. associaTR reads <STRn> as an allele of length n
                # (given RL=0 and a 1-char RU) so summed dosage = CN.
                for _n in range(1, 201):
                    fout.write(f'##ALT=<ID=STR{_n}>\n')
                # associaTR requires sample IDs to be strictly numeric, so we
                # remove the 'CPG' prefix
                fout.write(line.replace('CPG', ''))

            elif line.startswith('##'):
                # Pass meta-information lines through, but drop any INFO/FORMAT
                # header whose ID we redefine below (see _REDEFINED_HEADER_IDS)
                # so the output has exactly one definition per ID.
                if _is_redefined_header(line):
                    continue
                fout.write(line)

            else:
                # Process variant lines
                parts = line.strip().split('\t')
                n_seen += 1

                chrom = parts[0]
                if not chrom.startswith('chr'):
                    parts[0] = 'chr' + chrom

                record_id = parts[2] if parts[2] not in ('', '.') else 'NA'
                filter_field = parts[6]
                info_field = parts[7]
                format_field = parts[8]
                sample_data = parts[9:]

                is_multiallelic = 'MULTIALLELIC' in filter_field.split(';')

                info = parse_info_field(info_field)

                # Pull real END/SVTYPE from the GATK-SV INFO field rather than
                # faking them -- these are meaningful for SVs (unlike the SNP
                # version, where END is just set equal to POS).
                end = info.get('END', parts[1])
                svtype = info.get('SVTYPE', 'UNKNOWN')

                # RU stores ref/alt (here, the SV's symbolic alleles) so we can
                # discriminate between loci sharing the same CHR:POS
                new_info_fields = [
                    f'END={end}',
                    'REF=3',
                    'REPID=.',
                    'RL=0',
                    f'RU={parts[3]}-{parts[4]}',
                    f'SVTYPE={svtype}',
                    'VARID=.',
                ]
                updated_info_field = ';'.join(new_info_fields)

                # Update FORMAT field
                updated_format_field = 'GT:ADFL:ADIR:ADSP:LC:REPCI:REPCN:SO:QUAL'

                # Locate GT within the original FORMAT string by name, since
                # GATK-SV FORMAT has several subfields
                # (GT:ECN:EV:GQ:OGQ:RD_CN:RD_GQ:SL) and GT isn't guaranteed to
                # be the only field even though it is field 0 by VCF spec --
                # this makes the lookup explicit rather than assumed.
                format_keys = format_field.split(':')

                # Build the mock per-sample GT strings and the paired
                # per-sample dosage values used for the monomorphic drop.
                # MULTIALLELIC CNVs go through the CN path (emit "CN/0");
                # everything else goes through the standard GT path.
                mock_gts = []
                summed_gt_values = []

                if is_multiallelic:
                    if 'CN' not in format_keys:
                        # No CN available for a multi-allelic CNV -> no dosage
                        # to recover; drop the site.
                        drop_counts['multiallelic_no_cn'] += 1
                        _log_drop(
                            'multiallelic_no_cn',
                            parts[0], parts[1], record_id, svtype,
                            extra=f'FORMAT={format_field}',
                        )
                        continue
                    cn_index = format_keys.index('CN')
                    for sample in sample_data:
                        sample_parts = sample.split(':')
                        cn_val = sample_parts[cn_index] if cn_index < len(sample_parts) else '.'
                        if cn_val in ('', '.'):
                            mock_gts.append('./.')
                            continue
                        try:
                            cn_int = int(cn_val)
                        except ValueError:
                            mock_gts.append('./.')
                            continue
                        mock_gts.append(f'{cn_int}/0')
                        summed_gt_values.append(cn_int)
                else:
                    gt_index = format_keys.index('GT') if 'GT' in format_keys else 0
                    for sample in sample_data:
                        sample_parts = sample.split(':')
                        gt = sample_parts[gt_index].replace('|', '/')
                        mock_gts.append(gt if gt not in ('', '.') else './.')
                        if gt in ('', '.', './.'):
                            continue
                        alleles = gt.split('/')
                        if any(a == '.' for a in alleles):
                            continue
                        try:
                            summed_gt_values.append(sum(int(a) for a in alleles))
                        except ValueError:
                            continue

                if len(summed_gt_values) == 0:
                    drop_counts['monomorphic_no_callable'] += 1
                    _log_drop(
                        'monomorphic_no_callable',
                        parts[0], parts[1], record_id, svtype,
                        extra=f'n_samples={len(sample_data)}',
                    )
                    continue
                if all(sum_gt == summed_gt_values[0] for sum_gt in summed_gt_values):
                    drop_counts['monomorphic_uniform_dosage'] += 1
                    _log_drop(
                        'monomorphic_uniform_dosage',
                        parts[0], parts[1], record_id, svtype,
                        extra=f'dosage={summed_gt_values[0]};n_callable={len(summed_gt_values)}',
                    )
                    continue

                updated_sample_data = [':'.join([gt] + ['.'] * 8) for gt in mock_gts]
                n_written += 1

                # Update ALT column. MULTIALLELIC CN records emit GTs like "3/0",
                # "4/0", ... so the ALT column must declare enough symbolic
                # alleles for those indices to resolve; otherwise associaTR
                # crashes with IndexError inside GetLengthGenotypes (allele_lens
                # lookup out of bounds). associaTR reads <STRn> as an allele of
                # length n, so ALT = <STR1>,<STR2>,...,<STRmax_cn> gives:
                #   GT 0/0 -> summed dosage 0
                #   GT k/0 -> summed dosage k   (equal to CN)
                if is_multiallelic:
                    max_cn = max(summed_gt_values)  # >=1 here (0-only would have been dropped as monomorphic)
                    parts[4] = ','.join(f'<STR{i}>' for i in range(1, max_cn + 1))
                else:
                    parts[4] = '<STR1>'

                # Write the updated line to the output file
                updated_line = '\t'.join(
                    parts[:7] + [updated_info_field, updated_format_field] + updated_sample_data,
                )
                fout.write(updated_line + '\n')

    total_dropped = sum(drop_counts.values())
    print(
        f'SUMMARY\tinput={vcf_file_path}\tseen={n_seen}\twritten={n_written}\t'
        f'dropped={total_dropped}\t'
        + '\t'.join(f'{k}={v}' for k, v in drop_counts.items()),
        file=sys.stderr,
    )

    output_file.upload_from('temporary_gt_file.vcf.gz')


@click.option('--vcf-path', required=True, help='Input VCF file')
@click.option('--job-storage', default='20G')
@click.option('--job-cpu', default=1)
@click.command()
def main(vcf_path, job_storage, job_cpu):
    b = get_batch(name='SV VCF maker for associaTR')

    sv_vcf = vcf_path
    output_file = output_path(
        f'sv_vcf_for_associatr/{os.path.basename(vcf_path)}',
    )

    reformatting_job = b.new_python_job(name='Reformatting VCF')
    reformatting_job.storage(job_storage)
    reformatting_job.cpu(job_cpu)
    reformatting_job.call(reformat_vcf, sv_vcf, output_file)

    b.run(wait=False)


if __name__ == '__main__':
    main()
