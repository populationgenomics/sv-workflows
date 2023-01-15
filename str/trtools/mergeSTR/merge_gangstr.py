#!/usr/bin/env python

"""
Running Batch script to merge GangSTR vcf.gz files into one combined VCF. 
"""
import os
import hailtop.batch as hb
import click

DATASET = os.getenv('DATASET')
HAIL_BUCKET = os.getenv('HAIL_BUCKET')
OUTPUT_SUFFIX = os.getenv('OUTPUT')
BILLING_PROJECT = os.getenv('HAIL_BILLING_PROJECT')
ACCESS_LEVEL = os.getenv('ACCESS_LEVEL')

REF_FASTA = 'gs://cpg-reference/hg38/v1/Homo_sapiens_assembly38.fasta'
TRTOOLS_IMAGE = "australia-southeast1-docker.pkg.dev/cpg-common/images/trtools:v4.0.2"
BCFTOOLS_IMAGE = "australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.10.2--h4f4756c_2"


# cram = b.read_input_group(**{'cram': cram_path, 'cram.crai': cram_path + '.crai'})
@click.command()

def main():  # pylint: disable=missing-function-docstring 

    # Initializing Batch
    backend = hb.ServiceBackend(billing_project=BILLING_PROJECT, bucket=HAIL_BUCKET)
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    input_file = b.read_input_group(
    TOB01784 = "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB01784_GangSTR.vcf.gz",
    TOB01784_tbi = "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB01784_GangSTR.vcf.gz.tbi",
    TOB01791= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB01791_GangSTR.vcf.gz",
    TOB01791_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB01791_GangSTR.vcf.gz.tbi",
    TOB0901= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB0901_GangSTR.vcf.gz",
    TOB0901_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB0901_GangSTR.vcf.gz.tbi",
    TOB1086= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1086_GangSTR.vcf.gz",
    TOB1086_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1086_GangSTR.vcf.gz.tbi",
    TOB1170= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1170_GangSTR.vcf.gz",
    TOB1170_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1170_GangSTR.vcf.gz.tbi",
    TOB1213= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1213_GangSTR.vcf.gz",
    TOB1213_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1213_GangSTR.vcf.gz.tbi",
    TOB1458= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1458_GangSTR.vcf.gz",
    TOB1458_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1458_GangSTR.vcf.gz.tbi",
    TOB1567= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1567_GangSTR.vcf.gz",
    TOB1567_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1567_GangSTR.vcf.gz.tbi",
    TOB1578= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1578_GangSTR.vcf.gz",
    TOB1578_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1578_GangSTR.vcf.gz.tbi",
    TOB1751= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1751_GangSTR.vcf.gz",
    TOB1751_tbi= "gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/TOB1751_GangSTR.vcf.gz.tbi"
)

    trtools_job = b.new_job(name = "mergeSTR")
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('200G')
    trtools_job.cpu(16)

    trtools_job.declare_resource_group(ofile = {'vcf': '{root}.vcf'})
        
    trtools_job.command(f"""
     
    mergeSTR --vcfs {input_file.TOB01784},{input_file.TOB01791},{input_file.TOB0901},{input_file.TOB1086},{input_file.TOB1170},{input_file.TOB1213},{input_file.TOB1458},{input_file.TOB1567},{input_file.TOB1578},{input_file.TOB1751} --out {trtools_job.ofile} --vcftype gangstr
     
    """)

    gangstr_out_name = "tob10_GangSTR_merge"
    gangstr_output_path = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/merge/{gangstr_out_name}'
    b.write_output(trtools_job.ofile, gangstr_output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
