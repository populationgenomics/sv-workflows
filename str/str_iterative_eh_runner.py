#!/usr/bin/env python3

"""
Running Batch script to use ExpansionHunter to call STRs on 10 TOB-WGS genomes. 
"""
import os
import hailtop.batch as hb
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    remote_tmpdir,
)

config = get_config()

DATASET = config['workflow']['dataset']
HAIL_BUCKET = config['
OUTPUT_SUFFIX = config['workflow']['output_prefix']
BILLING_PROJECT = config['hail']['billing_project']
REF_FASTA = os.path.join(config['workflow']['reference_prefix'], 'hg38/v0/Homo_sapiens_assembly38.fasta')
SAMTOOLS_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/samtools:v0'
EH_IMAGE = "australia-southeast1-docker.pkg.dev/cpg-common/images/expansionhunter:5.0.0"

#Variant catalogs (required input)
EH_regions_path = 'gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json'

input_cram_dict={
    "TOB01784": "gs://cpg-tob-wgs-main/cram/nagim/CPG3160.cram",
    "TOB01791": "gs://cpg-tob-wgs-main/cram/nagim/CPG3210.cram",
    "TOB0901": "gs://cpg-tob-wgs-main/cram/nagim/CPG8458.cram",
    "TOB1086": "gs://cpg-tob-wgs-main/cram/nagim/CPG9688.cram",
    "TOB1170": "gs://cpg-tob-wgs-main/cram/nagim/CPG3954.cram",
    "TOB1213": "gs://cpg-tob-wgs-main/cram/nagim/CPG4358.cram",
    "TOB1458": "gs://cpg-tob-wgs-main/cram/nagim/CPG6551.cram",
    "TOB1567": "gs://cpg-tob-wgs-main/cram/nagim/CPG1503.cram",
    "TOB1578": "gs://cpg-tob-wgs-main/cram/nagim/CPG1610.cram",
    "TOB1751": "gs://cpg-tob-wgs-main/cram/nagim/CPG3053.cram"
}

@click.command()

def main():  # pylint: disable=missing-function-docstring 
   # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    EH_regions = b.read_input(EH_regions_path)

    #Iterate over each sample and perform 3 jobs 1) Index CRAM 2) Expansion Hunter
    for cram in list(input_cram_dict.keys()):

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(**{'cram': input_cram_dict[cram], 'cram.crai': input_cram_dict[cram]+ '.crai'})

        #Samtools job initialisation
        samtools_job = b.new_job(name = f'Index {cram}')
        samtools_job.image(SAMTOOLS_IMAGE)
        samtools_job.storage('200G')
        samtools_job.cpu(16)
        samtools_job.command(f"""samtools index -@ 20 {crams['cram']}""")

        # Working with CRAM files requires the reference fasta
        ref = b.read_input_group(
            **dict(
                base=REF_FASTA,
                fai=REF_FASTA + '.fai',
                dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
                + '.dict',
            )
        )
        
        # ExpansionHunter job initialisation
        eh_job = b.new_job(name = f'{cram}_EH')
        eh_job.image(EH_IMAGE)
        eh_job.depends_on(samtools_job)
        eh_job.storage('200G')
        eh_job.cpu(8)

        eh_job.declare_resource_group(ofile = {'vcf': '{root}.vcf',
                                               'json': '{root}.json',
                                               'realigned_bam': '{root}_realigned.bam'
        })

        eh_job.command(f"""
        ExpansionHunter --reads {crams['cram']} --reference {ref.base} --variant-catalog {EH_regions} --threads 16 --analysis-mode streaming --output-prefix {eh_job.ofile}
        """
        )
        # ExpansionHunter output writing
        eh_out_fname = f'{cram}_EH'
        eh_output_path = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/{eh_out_fname}'
        b.write_output(eh_job.ofile, eh_output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
