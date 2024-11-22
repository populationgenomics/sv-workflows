#!/usr/bin/env python3

"""
This script will

- extract a SNV VCFs from a VDS.

analysis-runner \
   --description "get common and rare variant VCFs" \
   --dataset "bioheart" \
   --access-level "test" \
   --output-dir raw-vds/hope-backup-vcfs \
    python3 get_genotype_vcf.py --vds-path=gs://cpg-bioheart-test/vds/tenk10k1-0.vds --chromosomes chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
"""

import logging

import click
import pandas as pd

import hail as hl
from hail.methods import export_vcf

from cpg_utils import to_path
from cpg_utils.config import get_config, output_path
from cpg_utils.hail_batch import get_batch, init_batch


def can_reuse(path: str):
    """
    checks for existence of a Path
    if the path is a MT or VDS, checks for the success file
    """

    if not path:
        return False

    path_as_path = to_path(path)

    if path_as_path.suffix in ['.mt', '.ht']:
        path_as_path /= '_SUCCESS'
    if path_as_path.suffix in ['.vds']:
        path_as_path /= 'variant_data/_SUCCESS'

    return path_as_path.exists()


def add_remove_chr_and_index_job(vcf_path):
    """
    Reads a VCF file, it creates an .csi index file and
    and removes "chr" from the chromosome values

    Args:
    vcf_path: input path
    """
    # remove chr & add index file using bcftools
    vcf_input = get_batch().read_input(vcf_path)

    vcf_size = to_path(vcf_path).stat().st_size
    storage_required = ((vcf_size // 1024**3) or 1) * 2.2
    bcftools_job = get_batch().new_job(name='remove chr and index vcf')
    bcftools_job.declare_resource_group(
        output={
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.csi': '{root}.vcf.bgz.csi',
        },
    )
    bcftools_job.image(get_config()['images']['bcftools'])
    bcftools_job.cpu(4)
    bcftools_job.storage(f'{storage_required}Gi')
    # now remove "chr" from chromosome names using bcftools
    bcftools_job.command('for num in {1..22} X Y; do echo "chr${num} ${num}" >> chr_update.txt; done')
    bcftools_job.command(
        f"""
        bcftools annotate --rename-chrs chr_update.txt --set-id +'%CHROM\:%POS\:%REF\:%FIRST_ALT' {vcf_input} | \\
        bgzip -c > {bcftools_job.output['vcf.bgz']}
        bcftools index -c {bcftools_job.output['vcf.bgz']}
    """,
    )
    logging.info('VCF rename/index jobs scheduled!')

    # save both output files
    get_batch().write_output(bcftools_job.output, vcf_path.removesuffix('.vcf.bgz'))


# inputs:
@click.option('--vds-path', help=' GCP gs:// path to the VDS')
@click.option('--chromosomes', help=' e.g., chr22,chrX ')
@click.option('--cv-maf-threshold', default=0.01)
@click.command()
def main(
    vds_path: str,
    chromosomes: str,
    cv_maf_threshold: float,
):
    """
    Write genotypes as VCF
    """

    init_batch(worker_memory='highmem')

    vds = hl.vds.read_vds(vds_path)
    vds_name = vds_path.split('/')[-1].split('.')[0]

    for chromosome in chromosomes.split(','):
        # create paths and check if they exist already
        cv_vcf_path = output_path(f'vds-{vds_name}/{chromosome}_common_variants.vcf.bgz')
        cv_vcf_existence_outcome = can_reuse(cv_vcf_path)
        logging.info(f'Does {cv_vcf_path} exist? {cv_vcf_existence_outcome}')

        if not cv_vcf_existence_outcome:
            # consider only relevant chromosome
            chrom_vds = hl.vds.filter_chromosomes(vds, keep=chromosome)

            # split multiallelic loci (necessary pre-densifying)
            chrom_vds = hl.vds.split_multi(chrom_vds, filter_changed_loci=True)

            # densify to matrix table object
            mt = hl.vds.to_dense_mt(chrom_vds)

            # filter out loci & variant QC
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # remove hom-ref

            mt = hl.variant_qc(mt)

            if not cv_vcf_existence_outcome:
                # common variants only
                cv_mt = mt.filter_rows(hl.min(mt.variant_qc.AF) >= cv_maf_threshold)

                if not cv_vcf_existence_outcome:
                    # remove fields not in the VCF
                    cv_mt = cv_mt.drop('gvcf_info')

                    # export to vcf common variants only
                    export_vcf(cv_mt, cv_vcf_path)

        # check existence of index file (CV) separately
        cv_index_file_existence_outcome = can_reuse(f'{cv_vcf_path}.csi')
        logging.info(f'Does {cv_vcf_path}.csi exist? {cv_index_file_existence_outcome}')
        if not cv_index_file_existence_outcome:
            add_remove_chr_and_index_job(cv_vcf_path)

    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()  # pylint: disable=no-value-for-parameter
