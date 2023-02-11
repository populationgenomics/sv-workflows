#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses HipSTR to call STRs on WGS cram files.
For example:
analysis-runner --access-level test --dataset hgdp --description 'hipstr run' --output-dir 'str/sensitivity-analysis/hipstr/untrimmed_coordinates_1_based' hipstr_runner_hgdp.py --variant-catalog=gs://cpg-hgdp-test/str/untrimmed_coordinates_resources/hg38_hipstr_catalog_untrimmed_coordinates.bed --input-dir=gs://cpg-hgdp-test/cram/nagim

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import logging
import hailtop.batch as hb


import click
from sample_metadata.apis import SampleApi, ParticipantApi
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, reference_path
from cpg_workflows.batch import get_batch
from google.cloud import storage

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
HIPSTR_IMAGE = config['images']['hipstr']
BCFTOOLS_IMAGE = config['images']['bcftools']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
@click.option('--input-dir', help='Full path to directory containing crams')
@click.command()
def main(
    variant_catalog, input_dir):  # pylint: disable=missing-function-docstring
    # Initializing Batch 
    b = get_batch()

    ref_fasta = ('gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta')
  
    #ref_fasta = str(reference_path('broad/ref_fasta'))
    crams_gcp_path =[]
    data = ParticipantApi().get_participants('hgdp-test')
    population_groups = ["French", "Yoruba", "Han", "Northern Han"]
    participant_ids = []

    for i in data: 
        if i["meta"]["Population name"] in population_groups: 
            participant_ids.append(i['id']) 
    samples = SampleApi().get_samples(body_get_samples={'project_ids':['hgdp-test'], "participant_ids": participant_ids})
    for i in samples: 
        crams_gcp_path.append("gs://cpg-hgdp-test/cram/nagim/"+i["id"]+".cram")
    
    hipstr_regions = b.read_input(variant_catalog)

    cram_collection = {}

    # Iterate over each sample to add to cram_collection
    for cram_obj in crams_gcp_path:
        cpg_sample_id = cram_obj.replace('.cram','')[30:]
        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj, 'cram.crai': cram_obj + '.crai'}
        )
        cram_collection[cpg_sample_id] = crams['cram']
    
    crams_batch_path = ""
    for i in cram_collection:
        crams_batch_path+=str(cram_collection[i])
        crams_batch_path+= ","
    crams_batch_path = crams_batch_path[:-1]

    ref = b.read_input_group(
            **dict(
                base=ref_fasta,
                fai=ref_fasta + '.fai',
                dict=ref_fasta.replace('.fasta', '')
                .replace('.fna', '')
                .replace('.fa', '')
                + '.dict',
            )
        )

    hipstr_job = b.new_job(name=f'HipSTR running')
    hipstr_job.image(HIPSTR_IMAGE)
    hipstr_job.storage('1500G')
    hipstr_job.cpu(16)

    hipstr_job.declare_resource_group(
        hipstr_output={
            'vcf.gz': '{root}.vcf.gz',
            'viz.gz': '{root}.viz.gz', 
            'log.txt': '{root}.log.txt'
        }
    )

    hipstr_job.command(
        f"""
    HipSTR --bams {crams_batch_path} --fasta {ref.base} --regions {hipstr_regions} --str-vcf {hipstr_job.hipstr_output['vcf.gz']} --viz-out {hipstr_job.hipstr_output['viz.gz']} --log {hipstr_job.hipstr_output['log.txt']} --output-filters
    """
    )
    # HipSTR output writing
    file_name = "91_hgdp_genomes"
    hipstr_output_path_vcf = output_path(f'{file_name}.vcf.gz')
    b.write_output(hipstr_job.hipstr_output['vcf.gz'], hipstr_output_path_vcf)

    hipstr_output_path_viz = output_path(f'{file_name}_hipstr.viz.gz')
    b.write_output(hipstr_job.hipstr_output['viz.gz'], hipstr_output_path_viz)

    hipstr_output_path_log = output_path(f'{file_name}_hipstr.log.txt')
    b.write_output(hipstr_job.hipstr_output['log.txt'], hipstr_output_path_log)

    samtools_job = b.new_job(name=f'tester Unzip VCF')
    samtools_job.image(SAMTOOLS_IMAGE)
    samtools_job.storage('20G')
    samtools_job.cpu(8)
    samtools_job.declare_resource_group(
        vcf = {
            'vcf': '{root}.vcf'
        }
    )

    samtools_job.command(
            f"""
            bgzip -d -c {hipstr_job.hipstr_output['vcf.gz']} > {samtools_job.vcf['vcf']}
        
            """
        )
    samtools_job_output_path = output_path(f'{file_name}_hipstr.vcf')
    b.write_output(samtools_job.vcf['vcf'], samtools_job_output_path)


    b.run(wait = False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
