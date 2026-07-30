#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
Cluster SVs across two GATK-SV cohort VCFs (BioHeart + TOB) and extract the
set of variants with carriers in BOTH cohorts.

Two-job pipeline:
  1. gatk SVCluster on the two cohort VCFs -> joint_clustered.vcf.gz
  2. bcftools split-by-cohort + isec -> common_sv.vcf.gz

Inputs (all in GCS):
  --bioheart-vcf       BioHeart cohort VCF (bgzipped + tabix-indexed)
  --tob-vcf            TOB cohort VCF (bgzipped + tabix-indexed)
  --bioheart-ploidy    BioHeart ploidy table (TSV, GATK-SV format)
  --tob-ploidy         TOB ploidy table (TSV, GATK-SV format)
  --reference-fasta    GRCh38 fasta (expects .fai and .dict siblings)

analysis-runner --dataset tenk10k-sv --access-level test \
    --output-dir pub-analysis/sv-cluster-common \
    --description "SVCluster across BioHeart + TOB" \
    sv_cluster_common.py \
        --bioheart-vcf=gs://cpg-tenk10k-test/pub-analysis/final-vcf/filtered/bioheart_common_maf_gte_1pct.vcf.gz \
        --tob-vcf=gs://cpg-tenk10k-test/pub-analysis/final-vcf/filtered/tob_common_maf_gte_1pct.vcf.gz \
        --ploidy=gs://cpg-tenk10k-test/sv/SVCluster/combined_ploidy.tsv \
        --reference-fasta=gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta
"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path


@click.command()
@click.option('--bioheart-vcf', required=True, help='BioHeart cohort VCF (bgzipped, .tbi sibling required)')
@click.option('--tob-vcf', required=True, help='TOB cohort VCF (bgzipped, .tbi sibling required)')
@click.option('--ploidy', required=True, help='Ploidy table (TSV, GATK-SV format)')
@click.option(
    '--reference-fasta',
    required=True,
    help='GRCh38 fasta; .fai and .dict siblings expected',
)
@click.option('--cluster-storage', default='100G', help='Disk for the SVCluster job')
@click.option('--cluster-memory', default='32G', help='Memory for the SVCluster job')
@click.option('--cluster-cpu', default=4, type=int, help='CPUs for the SVCluster job')
@click.option(
    '--cluster-algorithm',
    default='SINGLE_LINKAGE',
    type=click.Choice(['SINGLE_LINKAGE', 'MAX_CLIQUE', 'DEFRAGMENT_CNV']),
    help='SVCluster --algorithm',
)
def main(
    bioheart_vcf,
    tob_vcf,
    ploidy,
    reference_fasta,
    cluster_storage,
    cluster_memory,
    cluster_cpu,
    cluster_algorithm,
):
    b = get_batch(name='SVCluster BioHeart + TOB -> common SVs')

    # Read inputs. VCFs and reference come in as groups so the sibling indexes
    # land next to the primary file inside the job.
    bh_vcf = b.read_input_group(base=bioheart_vcf, tbi=bioheart_vcf + '.tbi')
    tb_vcf = b.read_input_group(base=tob_vcf, tbi=tob_vcf + '.tbi')
    r_ploidy = b.read_input(ploidy)
    ref = b.read_input_group(
        fasta=reference_fasta,
        fai=reference_fasta + '.fai',
        dict=reference_fasta.replace('.fasta', '.dict'),
    )

    # ---- Job 1: SVCluster --------------------------------------------------
    cluster_job = b.new_job(name='SVCluster BioHeart + TOB')
    cluster_job.image(get_config()['images']['gatk'])
    cluster_job.storage(cluster_storage)
    cluster_job.memory(cluster_memory)
    cluster_job.cpu(cluster_cpu)
    cluster_job.declare_resource_group(
        clustered={
            'vcf.gz': '{root}.vcf.gz',
            'vcf.gz.tbi': '{root}.vcf.gz.tbi',
        },
    )
    cluster_job.command(
        f"""
        set -euxo pipefail

        gatk --java-options "-Xmx{cluster_memory.rstrip('G')}g" SVCluster \\
            -V {bh_vcf.base} \\
            -V {tb_vcf.base} \\
            -O {cluster_job.clustered['vcf.gz']} \\
            -R {ref.fasta} \\
            --ploidy-table {r_ploidy} \\
            --algorithm {cluster_algorithm}
        """,
    )

    b.write_output(cluster_job.clustered, output_path('joint_clustered'))


    b.run(wait=False)


if __name__ == '__main__':
    main()
