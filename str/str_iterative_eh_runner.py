#!/usr/bin/env python3

"""
Running Batch script to use ExpansionHunter to call STRs on 10 TOB-WGS genomes. 
"""
import os
import hailtop.batch as hb
import click

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.apis import AnalysisApi

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

#required input

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

    tob_wgs_to_cpg_sample_id_map: dict[str, str] = SampleApi().get_sample_id_map_by_external(tob_wgs_ids, project=DATASET)
    cpg_sample_id_to_tob_wgs_id = {cpg_id: tob_wgs_id for tob_wgs_id, cpg_id in tob_wgs_to_cpg_sample_id_map.items()}

    analysis_query_model = AnalysisQueryModel(
        sample_ids=tob_wgs_to_cpg_sample_id_map.values(),
        projects=[DATASET],
        type=AnalysisType("cram"),
        status=AnalysisStatus("completed"),
        meta={"sequence_type": "genome", "source": "nagim"}
    )
    crams_path = AnalysisApi().query_analyses(analysis_query_model)

    for cram_object in crams:
        cram = cram_object["output"]
        cpg_sample_id = cram_object["sample_ids"][0]
        tob_wgs_id = cpg_sample_id_to_tob_wgs_id[cpg_sample_id]
    EH_regions = b.read_input(ehregions)

    #Iterate over each sample and perform 3 jobs 1) Index CRAM 2) Expansion Hunter
    for cram in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(**{'cram': cram, 'cram.crai': cram+ '.crai'})

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
        eh_output_path = output_path(f'{cram}_EH')
        b.write_output(eh_job.ofile, eh_output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
