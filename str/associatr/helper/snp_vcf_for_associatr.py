#!/usr/bin/env python3

"""
This script is used to convert existing SNP VCF files (chromosome-specific) into a mock ExpansionHunter-style VCF for direct use with associaTR.

Please note ad-hoc changes:
- RU field stores the REF and ALT alleles in the format of 'REF-ALT'
- RL field is set to 0 for all loci so that the REF allele is coded as 0 (https://github.com/gymrek-lab/TRTools/blob/master/trtools/utils/tr_harmonizer.py#L515)
- REF field is set to 3 (arbitrary value)

analysis-runner --dataset bioheart --access-level test --output-dir str/associatr --description "snp vcf for associatr" \
snp_vcf_for_associatr.py --vcf-dir=gs://cpg-bioheart-main/saige-qtl/bioheart_n990_and_tob_n1055/input_files/genotypes/vds-tenk10k1-0 \
--chromosomes=11,12,13,14,16,6,9 --job-storage=100G --job-cpu=8

"""

import gzip

import click

from cpg_utils import to_path
from cpg_utils.config import output_path
from cpg_utils.hail_batch import get_batch


def reformat_vcf(vcf_file_path, output_file_path):
    vcf_file = to_path(vcf_file_path)
    output_file = to_path(output_file_path)
    with gzip.open(vcf_file, 'rt') as fin, open('temporary_gt_file.txt', 'w') as fout:
        for line in fin:
            if line.startswith('##hailversion'):
                # Write header line and update FORMAT column
                fout.write(line)
                # Add new FORMAT fields to the header (ExpansionHunter style)
                fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n')
                fout.write('##FORMAT=<ID=ADFL,Number=1,Type=String,Description="">\n')
                fout.write('##FORMAT=<ID=ADIR,Number=1,Type=String,Description="">\n')
                fout.write('##FORMAT=<ID=ADSP,Number=1,Type=String,Description="">\n')
                fout.write('##FORMAT=<ID=LC,Number=1,Type=Float,Description="">\n')
                fout.write('##FORMAT=<ID=REPCI,Number=1,Type=String,Description="">\n')
                fout.write('##FORMAT=<ID=REPCN,Number=1,Type=String,Description="">\n')
                fout.write('##FORMAT=<ID=SO,Number=1,Type=String,Description="">\n')
                fout.write('##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="">\n')
                fout.write('##INFO=<ID=END,Number=1,Type=Integer,Description="">\n')
                fout.write('##INFO=<ID=REF,Number=1,Type=Integer,Description="">\n')
                fout.write('##INFO=<ID=REPID,Number=1,Type=String,Description="">\n')
                fout.write('##INFO=<ID=RL,Number=1,Type=Integer,Description="">\n')
                fout.write(
                    '##INFO=<ID=RU,Number=1,Type=String,Description="Storing the REF/ALT info instead of Repeat Unit to retain the REF/ALT info in the association output files, which is crucial where there are multiple variants with the same CHR:POS coordinates (eg multi allelic loci)">\n',
                )
                fout.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">\n')
                fout.write('##INFO=<ID=VARID,Number=1,Type=String,Description="">\n')
                fout.write('##ALT=<ID=STR1>\n')
            elif line.startswith('##'):
                # Write meta-information lines as they are
                fout.write(line)
            elif line.startswith('#CHROM'):
                # associaTR requires sample IDS to be strictly numeric so we remove 'CPG' prefix
                fout.write(line.replace('CPG', ''))

            else:
                # Process variant lines
                parts = line.strip().split('\t')
                parts[0] = 'chr' + parts[0]
                info_field = parts[7]
                format_field = parts[8]
                sample_data = parts[9:]

                # Extract PID from the sample data if available
                # pid_index = format_field.split(':').index('PID') if 'PID' in format_field.split(':') else None
                # pid = sample_data[0].split(':')[pid_index] if pid_index is not None else '.'

                # Update INFO field (add new fields required by ExpansionHunter, filler fields eg REF =0)
                # RU stores the ref and alt_allele; this will allow us to discriminate between loci with the same POS
                new_info_fields = [
                    f'END={parts[1]}',
                    'REF=3',
                    'REPID=.',
                    'RL=0',
                    f'RU={parts[3]}-{parts[4]}',
                    'VARID=.',
                ]
                updated_info_field = ';'.join(new_info_fields)

                # Update FORMAT field
                updated_format_field = 'GT:ADFL:ADIR:ADSP:LC:REPCI:REPCN:SO:QUAL'

                # Check if all GT sums are the same across all samples, excluding missing GT values
                summed_gt_values = []
                for sample in sample_data:
                    gt = sample.split(':')[0].replace('|', '/')
                    if gt != './.':
                        alleles = gt.split('/')
                        if alleles[0] != '.' and alleles[1] != '.':
                            gt_sum = int(alleles[0]) + int(alleles[1])
                            summed_gt_values.append(gt_sum)

                # Skip writing this line if all summed GT values are the same or if all GT values are missing
                if len(summed_gt_values) == 0 or all(sum_gt == summed_gt_values[0] for sum_gt in summed_gt_values):
                    continue

                # Update sample data
                updated_sample_data = []
                for sample in sample_data:
                    sample_parts = sample.split(':')
                    gt = sample_parts[0].replace('|', '/')  # Replace '|' with '/' in the GT value
                    new_sample = [gt] + ['.'] * 8
                    updated_sample_data.append(':'.join(new_sample))

                # Update ALT column
                parts[4] = '<STR1>'

                # Write the updated line to the output file
                updated_line = '\t'.join(parts[:7] + [updated_info_field, updated_format_field] + updated_sample_data)
                fout.write(updated_line + '\n')
    output_file.upload_from('temporary_gt_file.txt')


@click.option('--vcf-dir', required=True, help='Input VCF file')
@click.option(
    '--chromosomes',
    required=True,
    default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22',
    help='Chromosomes to process',
)
@click.option('--job-storage', default='20G')
@click.option('--job-cpu', default=1)
@click.command()
def main(vcf_dir, chromosomes, job_storage, job_cpu):
    b = get_batch(name='SNP VCF maker for associaTR')
    for chrom in chromosomes.split(','):
        snp_vcf = f'{vcf_dir}/chr{chrom}_common_variants.vcf.bgz'
        output_file = output_path(f'common_variants_snps/hail_filtered_chr{chrom}.vcf')
        reformatting_job = b.new_python_job(name=f'Reformating chr{chrom} VCF')
        reformatting_job.storage(job_storage)
        reformatting_job.cpu(job_cpu)
        reformatting_job.call(reformat_vcf, snp_vcf, output_file)
    b.run(wait=False)


if __name__ == '__main__':
    main()
