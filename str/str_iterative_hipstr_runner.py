#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses HipSTR to call STRs on WGS cram files, using the joint calling option. 
For example:
analysis-runner --access-level test --dataset hgdp --description 'hipstr run' --output-dir 'str/sensitivity-analysis/hipstr' str_iterative_hipstr_runner.py --variant-catalog=gs://cpg-hgdp-test/str/untrimmed_coordinates_resources/hg38_hipstr_catalog_untrimmed_coordinates.bed --output-file-name=hipster_90_genomes --dataset=hgdp HGDP00511

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import logging

import os
import click
from sample_metadata.apis import AnalysisApi, SampleApi
from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.models import AnalysisStatus
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
HIPSTR_IMAGE = config['images']['hipstr']
BCFTOOLS_IMAGE = config['images']['bcftools']


# inputs:
@click.option('--variant-catalog', help='Full path to HipSTR Variants catalog')
@click.option('--dataset', help='dataset eg tob-wgs')
@click.argument('external-wgs-ids', nargs=-1)
@click.option('--output-file-name', help='Output file name without file extension')
@click.command()
def get_cloudfuse_paths(dataset, external_wgs_ids):
    """Retrieves cloud fuse paths and outputs as a comma-separated string for crams associated with the external-wgs-ids"""
    # Extract CPG IDs from external sample IDs
    external_id_to_cpg_id: dict[str, str] = SampleApi().get_sample_id_map_by_external(
        dataset, list(external_wgs_ids)
    )
    analysis_query_model = AnalysisQueryModel(
        sample_ids=list(external_id_to_cpg_id.values()),
        projects=[dataset],
        type=AnalysisType('cram'),
        status=AnalysisStatus('completed'),
        meta={'sequence_type': 'genome', 'source': 'nagim'},
    )
    crams_path = AnalysisApi().query_analyses(analysis_query_model)
    cpg_sids_with_crams = set(sid for sids in crams_path for sid in sids['sample_ids'])
    cpg_id_to_external_id = {
        cpg_id: external_wgs_id
        for external_wgs_id, cpg_id in external_id_to_cpg_id.items()
    }
    cpg_sids_without_crams = set(cpg_id_to_external_id.keys()) - cpg_sids_with_crams
    if cpg_sids_without_crams:
        external_wgs_sids_without_crams = ', '.join(
            cpg_id_to_external_id[sid] for sid in cpg_sids_without_crams
        )
        logging.warning(
            f'There were some samples without CRAMs: {external_wgs_sids_without_crams}'
        )
    # Create string containing paths based on /cramfuse
    cramfuse_path = ''
    for cram_obj in crams_path:
        cpg_sample_id = os.path.basename(cram_obj['output']).split('.')[0]
        if dataset == 'tob-wgs':
            cramfuse_path += '/cramfuse/cram/nagim/' + cpg_sample_id + '.cram,'
        else:
            cramfuse_path += '/cramfuse/cram/' + cpg_sample_id + '.cram,'
    cramfuse_path = crams_path[:-1]  # removes the terminating comma
    return cramfuse_path


def main(
    variant_catalog, dataset, external_wgs_ids, output_file_name
):  # pylint: disable=missing-function-docstring
    b = get_batch()
    # Create HipSTR job
    hipstr_job = b.new_job(name=f'HipSTR running')
    hipstr_job.image(HIPSTR_IMAGE)
    hipstr_job.storage('375G')
    hipstr_job.cpu(16)

    hipstr_job.declare_resource_group(
        hipstr_output={
            'vcf.gz': '{root}.vcf.gz',
            'viz.gz': '{root}.viz.gz',
            'log.txt': '{root}.log.txt',
        }
    )
    hipstr_job.cloudfuse(f'cpg-{dataset}-main', '/cramfuse')
    ref_fasta = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'

    cramfuse_path = get_cloudfuse_paths(dataset, external_wgs_ids)

    # Read in HipSTR variant catalog
    hipstr_regions = b.read_input(variant_catalog)

    # Read in reference
    ref = b.read_input_group(
        **dict(
            base=ref_fasta,
            fai=ref_fasta + '.fai',
            dict=ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )

    hipstr_job.command(
        f"""
    HipSTR --bams {cramfuse_path} \\
        --fasta {ref.base} \\
        --regions {hipstr_regions} \\
        --str-vcf {hipstr_job.hipstr_output['vcf.gz']} \\
        --viz-out {hipstr_job.hipstr_output['viz.gz']} \\
        --log {hipstr_job.hipstr_output['log.txt']} \\
        --output-filters
    """
    )
    # HipSTR output writing
    hipstr_output_path_vcf = output_path(f'{output_file_name}.vcf.gz')
    b.write_output(hipstr_job.hipstr_output['vcf.gz'], hipstr_output_path_vcf)

    hipstr_output_path_viz = output_path(f'{output_file_name}_hipstr.viz.gz')
    b.write_output(hipstr_job.hipstr_output['viz.gz'], hipstr_output_path_viz)

    hipstr_output_path_log = output_path(f'{output_file_name}_hipstr.log.txt')
    b.write_output(hipstr_job.hipstr_output['log.txt'], hipstr_output_path_log)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
