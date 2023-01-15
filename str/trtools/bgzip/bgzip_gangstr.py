#!/usr/bin/env python

"""
Running Batch script to sort, bgzip, and index the GangSTR vcfs. 
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


input_vcf_dict_gangstr={
    "TOB01784": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB01784_GangSTR.vcf",
    "TOB01791": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB01791_GangSTR.vcf",
    "TOB0901": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB0901_GangSTR.vcf",
    "TOB1086": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1086_GangSTR.vcf",
    "TOB1170": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1170_GangSTR.vcf",
    "TOB1213": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1213_GangSTR.vcf",
    "TOB1458": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1458_GangSTR.vcf",
    "TOB1567": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1567_GangSTR.vcf",
    "TOB1578": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1578_GangSTR.vcf",
    "TOB1751": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1751_GangSTR.vcf"
}

@click.command()

def main():  # pylint: disable=missing-function-docstring 

    # Initializing Batch
    backend = hb.ServiceBackend(billing_project=BILLING_PROJECT, bucket=HAIL_BUCKET)
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    #Iterate over each sample and perform 3 jobs 1) Index CRAM 2) GangSTR 3) Expansion Hunter
    for vcf in list(input_vcf_dict_gangstr.keys()):
        
        bcftools_job = b.new_job(name = f'{vcf}Files prep')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.storage('200G')
        bcftools_job.cpu(16)

        bcftools_job.declare_resource_group(
        gang_vcf={'vcf.gz': '{root}.gang.vcf.gz', 'vcf.gz.tbi': '{root}.gang.vcf.gz.tbi'}
    )

        vcfs = b.read_input_group(**{'vcf': input_vcf_dict_gangstr[vcf]})

        bcftools_job.command(f"""

        bcftools sort {vcfs['vcf']} | bgzip -c > {bcftools_job.gang_vcf['vcf.gz']}

        tabix -f -p vcf {bcftools_job.gang_vcf['vcf.gz']} > {bcftools_job.gang_vcf['vcf.gz.tbi']}
        
        """)

        # GangSTR output writing 
        gangstr_out_fname_gz = f'{vcf}_GangSTR.vcf.gz'
        gangstr_output_path_gz = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/{gangstr_out_fname_gz}'
        b.write_output(bcftools_job.gang_vcf['vcf.gz'], gangstr_output_path_gz)
       
        gangstr_out_fname_tbi = f'{vcf}_GangSTR.vcf.gz.tbi'
        gangstr_output_path_tbi = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/gangstr/{gangstr_out_fname_tbi}'
        b.write_output(bcftools_job.gang_vcf['vcf.gz.tbi'], gangstr_output_path_tbi)


    b.run(wait=False)

if __name__ == '__main__':
    main() 
