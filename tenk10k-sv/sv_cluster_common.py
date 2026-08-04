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

    # ---- Preprocess: add missing FORMAT header declarations so htsjdk will
    # emit records that carry CNV-related per-sample fields. Only lines not
    # already present in the header are appended, so this is a no-op if the
    # source VCF is already complete.
    def reheader_job(name, vcf_group):
        j = b.new_job(name=f'Reheader {name}')
        j.image(get_config()['images']['bcftools'])
        j.storage('20G')
        j.cpu(2)
        j.declare_resource_group(
            fixed={
                'vcf.gz': '{root}.vcf.gz',
                'vcf.gz.tbi': '{root}.vcf.gz.tbi',
            },
        )
        j.command(
            f"""
            set -euxo pipefail

            bcftools view -h {vcf_group.base} > hdr.txt

            add_format() {{
                local id="$1" number="$2" type="$3" desc="$4"
                grep -q "^##FORMAT=<ID=${{id}}," hdr.txt || \\
                    sed -i "/^#CHROM/i ##FORMAT=<ID=${{id}},Number=${{number}},Type=${{type}},Description=\\"${{desc}}\\">" hdr.txt
            }}

            add_format CN   1 Integer "Copy number"
            add_format CNQ  1 Integer "Copy number genotype quality"
            add_format ECN  1 Integer "Expected copy number"
            add_format RD_CN 1 Integer "Read-depth-based copy number"
            add_format RD_GQ 1 Integer "Read-depth-based genotype quality"
            add_format SL   1 Integer "Second-best likelihood"
            add_format OGQ  1 Integer "Original genotype quality"
            add_format EV   . String  "Classes of evidence supporting the call"

            bcftools reheader -h hdr.txt {vcf_group.base} -o {j.fixed['vcf.gz']}
            bcftools index -t -o {j.fixed['vcf.gz.tbi']} {j.fixed['vcf.gz']}
            echo "Records: $(bcftools view -H {j.fixed['vcf.gz']} | wc -l)"
            """,
        )
        return j

    bh_filter = reheader_job('bioheart', bh_vcf)
    tb_filter = reheader_job('tob', tb_vcf)

    # ---- Job 1: SVCluster --------------------------------------------------
    cluster_job = b.new_job(name='SVCluster BioHeart + TOB')
    cluster_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images/sv/gatk:2025-05-20-4.6.2.0-4-g1facd911e-NIGHTLY-SNAPSHOT')
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
            -V {bh_filter.fixed['vcf.gz']} \\
            -V {tb_filter.fixed['vcf.gz']} \\
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
