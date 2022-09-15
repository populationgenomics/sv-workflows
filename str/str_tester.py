#!/usr/bin/env python3

"""
This script uses ExpansionHunterv5 to call STRs on WGS cram files.
Required input: --ehregions (file path to variant catalog) and external sample IDs 
For example: 
analysis-runner --access-level test --dataset tob-wgs --description "tester" --output-dir "tester" str_iterative_eh_runner.py --ehregions=gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json TOB1XXXX TOB1XXXX

Required packages: sample-metadata, hail, click, os 
"""
import os
import hailtop.batch as hb
import click

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.apis import AnalysisApi, SampleApi
from sample_metadata.models import AnalysisStatus
from cpg_utils.config import get_config

from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    remote_tmpdir,
    output_path
)

config = get_config()

DATASET = config['workflow']['dataset']
OUTPUT_SUFFIX = config['workflow']['output_prefix']
BILaLING_PROJECT = config['hail']['billing_project']
REF_FASTA = os.path.join(config['workflow']['reference_prefix'], 'hg38/v0/Homo_sapiens_assembly38.fasta')
SAMTOOLS_IMAGE = os.path.join(config['workflow']['image_registry_prefix'], 'samtools:v0')
EH_IMAGE = os.path.join(config['workflow']['image_registry_prefix'], 'expansionhunter:5.0.0')

#inputs: 
#variant catalog
@click.option('--ehregions', help='Full path to Illumina Variants catalog')\
#input TOB ID 
@click.argument('tob-wgs-ids', nargs=-1)
@click.command()

def main(ehregions, tob_wgs_ids: list[str]):  # pylint: disable=missing-function-docstring 
   # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    crams_path = ["gs://cpg-tob-wgs-test/cram/nagim/CPG199760.cram"]
    EH_regions = b.read_input(ehregions)

    #Iterate over each sample and perform 3 jobs 1) Index CRAM 2) Expansion Hunter
    for cram_obj in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(**{'cram': cram_obj, 'cram.crai': cram_obj+ '.crai'})
        cpg_sample_id = "CPG199760"

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
        eh_job = b.new_job(name = f'ExpansionHunter:{cpg_sample_id} running')
        eh_job.image(EH_IMAGE)
        eh_job.storage('30G')
        eh_job.cpu(8)

        eh_job.declare_resource_group(ofile = {'vcf': '{root}.vcf',
                                               'json': '{root}.json',
                                               'realigned_bam': '{root}_realigned.bam'
        })

        eh_job.command(f"""
        ExpansionHunter  \\
        --reads {crams['cram']} \\
        --reference {ref.base} --variant-catalog {EH_regions}\\
         --threads 16 --analysis-mode streaming \\
         --output-prefix {eh_job.eh_output}
        """
        )
        # ExpansionHunter output writing
        eh_output_path = output_path(f'{cpg_sample_id}_EH')
        b.write_output(eh_job.eh_output, eh_output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
