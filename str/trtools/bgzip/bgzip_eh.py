#!/usr/bin/env python

"""
Running Batch script to bgzip, and index the EH vcfs. 
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


input_vcf_dict_eh={
    "TOB01784": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB01784_EH.vcf",
    "TOB01791": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB01791_EH.vcf",
    "TOB0901": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB0901_EH.vcf",
    "TOB1086": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1086_EH.vcf",
    "TOB1170": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1170_EH.vcf",
    "TOB1213": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1213_EH.vcf",
    "TOB1458": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1458_EH.vcf",
    "TOB1567": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1567_EH.vcf",
    "TOB1578": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1578_EH.vcf",
    "TOB1751": "gs://cpg-tob-wgs-test/hoptan-str/tob10/TOB1751_EH.vcf"
}

@click.command()

def main():  # pylint: disable=missing-function-docstring 

    # Initializing Batch
    backend = hb.ServiceBackend(billing_project=BILLING_PROJECT, bucket=HAIL_BUCKET)
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    ref = b.read_input_group(
        **dict(
            base=REF_FASTA,
            fai=REF_FASTA + '.fai',
            dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )

    #Iterate over each sample and perform bgzip and indexing 
    for vcf in list(input_vcf_dict_eh.keys()):
        
        bcftools_job = b.new_job(name = f'{vcf} EH Files prep')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.storage('200G')
        bcftools_job.cpu(16)

        bcftools_job.declare_resource_group(
        eh_vcf={'vcf.gz': '{root}.eh.vcf.gz', 'reheader.vcf.gz': '{root}.eh.reheader.vcf.gz', 'vcf.gz.tbi': '{root}.eh.reheader.vcf.gz.tbi'}
    )

        vcfs = b.read_input_group(**{'vcf': input_vcf_dict_eh[vcf]})

        bcftools_job.command(f"""

        bgzip -c {vcfs['vcf']} > {bcftools_job.eh_vcf['vcf.gz']}
        
        bcftools reheader -f {ref.fai} -o {bcftools_job.eh_vcf['reheader.vcf.gz']} {bcftools_job.eh_vcf['vcf.gz']} 

        tabix -f -p vcf {bcftools_job.eh_vcf['reheader.vcf.gz']} > {bcftools_job.eh_vcf['vcf.gz.tbi']}
        
        """)

        # EH output writing 
        eh_out_fname_gz = f'{vcf}_EHreheader.vcf.gz'
        eh_output_path_gz = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/{eh_out_fname_gz}'
        b.write_output(bcftools_job.eh_vcf['reheader.vcf.gz'], eh_output_path_gz)
       
        eh_out_fname_tbi = f'{vcf}_EHreheader.vcf.gz.tbi'
        eh_output_path_tbi = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/bgzip/eh/{eh_out_fname_tbi}'
        b.write_output(bcftools_job.eh_vcf['vcf.gz.tbi'], eh_output_path_tbi)


    b.run(wait=False)

if __name__ == '__main__':
    main() 
