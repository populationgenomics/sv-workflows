#!/usr/bin/env python3

"""
This script uses ExpansionHunterv5 to call STRs on WGS cram files.
Required input: --ehregions (file path to variant catalog) and external sample IDs 
For example: 
analysis-runner --access-level test --dataset tob-wgs --description "tester" --output-dir "tester" str_iterative_eh_runner.py --ehregions=gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json TOB1XXXX TOB1XXXX

Required packages: sample-metadata, hail, click, os 
pip install sample-metadata hail click

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
BILLING_PROJECT = config['hail']['billing_project']
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

    tob_wgs_to_cpg_sample_id_map: dict[str, str] = SampleApi().get_sample_id_map_by_external('tob-wgs',list(tob_wgs_ids))
    cpg_sample_id_to_tob_wgs_id = {cpg_id: tob_wgs_id for tob_wgs_id, cpg_id in tob_wgs_to_cpg_sample_id_map.items()}

    analysis_query_model = AnalysisQueryModel(
        sample_ids=list(tob_wgs_to_cpg_sample_id_map.values()),
        projects=['tob-wgs'],
        type=AnalysisType("cram"),
        status=AnalysisStatus("completed"),
        meta={"sequence_type": "genome", "source": "nagim"}
    )
    crams_path = AnalysisApi().query_analyses(analysis_query_model)
    EH_regions = b.read_input(ehregions)

    #Iterate over each sample and perform 3 jobs 1) Index CRAM 2) Expansion Hunter
    for cram_obj in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(**{'cram': cram_obj["output"], 'cram.crai': cram_obj["output"]+ '.crai'})
        cpg_sample_id = cram_obj["sample_ids"][0]
        tob_wgs_id = cpg_sample_id_to_tob_wgs_id[cpg_sample_id]
        
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
