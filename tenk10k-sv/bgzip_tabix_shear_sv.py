#!/usr/bin/env python3

"""
Shears a whole-genome EH-style SV VCF (the output of
tenk10k-sv/sv_vcf_for_associatr.py) into per-chromosome bgzipped + tabixed VCFs
suitable for tenk10k-sv/associatr_runner_sv.py, which expects one
hail_filtered_{chromosome}.vcf.bgz + .tbi pair per chromosome.

Differs from str/helper/bgzip_tabix.py in that the input here is a single
whole-genome VCF (not one file per chromosome), so this helper does the
per-chromosome split as well as the bgzip + tabix. The split is done with
awk directly on the gzip stream, so the input does not need to be indexed
beforehand.

analysis-runner --access-level test --dataset tenk10k-sv \
    --description 'shear SV VCF into per-chrom bgzip+tabix' \
    --output-dir 'pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/tob_common_chr_specific' \
    bgzip_tabix_shear_sv.py \
    --vcf-path=gs://cpg-tenk10k-sv-test/pub-analysis/final-vcf/filtered/sv_vcf_for_associatr/tob_common_maf_gte_1pct.vcf.gz \
    --chromosomes=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


@click.command()
@click.option('--vcf-path', required=True, help='Whole-genome EH-style SV VCF (.vcf.gz)')
@click.option(
    '--chromosomes',
    default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22',
    help='Comma-separated chromosome names (without chr prefix)',
)
@click.option('--job-memory', default='4G', help='Job memory')
@click.option('--job-storage', default='20G', help='Job storage')
@click.option('--job-cpu', default=1, help='Job CPU')
def main(vcf_path, chromosomes, job_memory, job_storage, job_cpu):
    """
    Splits a whole-genome VCF into per-chromosome shards and writes each shard
    out as a bgzip+tabix pair to the analysis-runner output directory.
    """
    b = get_batch(name='Shear SV VCF by chromosome')
    vcf_input = b.read_input(vcf_path)

    for chrom in chromosomes.split(','):
        chrom_name = f'chr{chrom}'

        shear_job = b.new_job(name=f'{chrom_name} SHEAR + BGZIP + TABIX')
        shear_job.image(BCFTOOLS_IMAGE)
        shear_job.memory(job_memory)
        shear_job.storage(job_storage)
        shear_job.cpu(job_cpu)

        shear_job.declare_resource_group(
            vcf_output={
                'vcf.bgz': '{root}.vcf.bgz',
                'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
            },
        )

        # zcat + awk keeps every header line plus records whose CHROM equals
        # chr{N}; bgzip -c re-encodes as bgzip (not plain gzip) so tabix can
        # build an index over the block offsets.
        shear_job.command(
            f"""

            zcat {vcf_input} \
                | awk -v c={chrom_name} 'BEGIN{{FS="\\t"; OFS="\\t"}} /^#/ {{print; next}} $1==c {{print}}' \
                | bgzip -c > {shear_job.vcf_output['vcf.bgz']}
            tabix -f -p vcf {shear_job.vcf_output['vcf.bgz']}

            """,
        )
        b.write_output(shear_job.vcf_output, output_path(f'hail_filtered_{chrom_name}'))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
