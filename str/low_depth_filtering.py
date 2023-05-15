#!/usr/bin/env python3
# pylint: disable=import-error, too-many-locals

"""
This script removes any entries in individual sample ExpansionHunter VCFs with FILTER set to "LowDepth"
"""
import click
from cloudpathlib import GSPath
from cyvcf2 import VCFReader, Writer
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

REF_FASTA = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
BCFTOOLS_IMAGE = config['images']['bcftools']

config = get_config()

ref = b.read_input_group(
        **dict(
            base=REF_FASTA,
            fai=REF_FASTA + '.fai',
            dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )

def low_depth_filter(vcf_input):
    """
    Creates a VCF with "LowDepth" calls filtered out
    """

    if isinstance(file, GSPath):
            file.copy(file.name)
            file = file.name
    reader = VCFReader(file)
    fname = "CPG10017_filtered.vcf"
    w = Writer(fname, vcf)
    for variant in reader: 
        if variant.FILTER != "LowDepth":
            w.write_record(variant)
    w.close()
    vcf.close()
    return w

@click.command()
@click.option('--input-dir', help='Input directory for ExpansionHunter VCFs')
def main(input_dir):
    # Initializing Batch

    files = to_path(input_dir).glob('*.reheader.vcf.gz')

    for file in files: 
        



if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter