#!/usr/bin/env python3


"""
This script uses HipSTR to call STRs on WGS cram files, using the joint calling option, and a sharded catalog.
For example:
analysis-runner --access-level test --dataset tob-wgs --description 'hipstr run' --output-dir 'str/sensitivity-analysis/hipstr' str_iterative_hipstr_runner.py  --output-file-name=hipster_90_genomes  --variant-catalog=gs://.... --dataset=hgdp HGDP00511

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import logging

import click

from metamist.graphql import gql, query
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch


config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
HIPSTR_IMAGE = config['images']['hipstr']
BCFTOOLS_IMAGE = config['images']['bcftools']


def get_cloudfuse_paths(dataset, input_cpg_sids):
    """Retrieves cloud fuse paths and outputs as a comma-separated string for crams associated with the external-wgs-ids"""

    cram_retrieval_query = gql(
        """
        query MyQuery($dataset: String!,$input_cpg_sids: [String!]!) {
    project(name: $dataset) {
        sequencingGroups(id: {in_: $input_cpg_sids}) {
        id
        sample {
            externalId
        }
        analyses(type: {eq: "cram"}, active: {eq: true}) {
            output
            timestampCompleted
        }
        }
    }

    }
        """
    )
    response = query(
        cram_retrieval_query,
        variables={'dataset': dataset, 'input_cpg_sids': input_cpg_sids},
    )
    crams_by_id = {}
    for i in response['project']['sequencingGroups']:
        for cram in i['analyses']:
            # ignore archived CRAMs and ONT crams
            if cram['output'] is None:
                continue
            if 'archive' in cram['output']:
                continue
            if 'ont' in cram['output']:
                continue
            crams_by_id[i['id']] = cram

    if len(crams_by_id) != len(input_cpg_sids):
        cpg_sids_without_crams = set(input_cpg_sids) - set(crams_by_id.keys())
        logging.warning(
            f'There were some samples without CRAMs: {cpg_sids_without_crams}'
        )

    # Create string containing paths based on /cramfuse
    cramfuse_path = []
    for cram_obj in crams_by_id.values():
        suffix = cram_obj['output'].removeprefix('gs://').split('/', maxsplit=1)[1]
        cramfuse_path.append(f'/cramfuse/{suffix}')
    cramfuse_path = ','.join(cramfuse_path)  # string format for input into hipstr
    return cramfuse_path


# inputs:
@click.option(
    '--job-storage', help='Storage of the Hail batch job eg 30G', default='30G'
)
@click.option('--job-memory', help='Memory of the Hail batch job', default='highmem')
@click.option(
    '--variant-catalog',
    help='Full path to HipSTR Variants sharded catalog directory, and file prefix',
)
@click.option('--dataset', help='dataset eg tob-wgs')
@click.option('--num-shards', help='Number of catalog shards', default=1)
@click.argument('internal-cpg-ids', nargs=-1)
@click.option('--output-file-name', help='Output file name without file extension')
@click.command()
def main(
    job_storage,
    job_memory,
    num_shards,
    variant_catalog,
    dataset,
    internal_cpg_ids,
    output_file_name,
):  # pylint: disable=missing-function-docstring
    b = get_batch()
    ref_fasta = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
    # Read in reference
    ref = b.read_input_group(
        **dict(
            base=ref_fasta,
            fai=ref_fasta + '.fai',
            dict=ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )

    cramfuse_path = get_cloudfuse_paths(dataset, internal_cpg_ids)

    for i in range(int(num_shards)):
        # Create HipSTR job
        hipstr_job = b.new_job(name=f'HipSTR shard {i+1} running')
        hipstr_job.image(HIPSTR_IMAGE)
        hipstr_job.storage(job_storage)
        hipstr_job.cpu(4)
        hipstr_job.memory(job_memory)

        hipstr_job.declare_resource_group(
            hipstr_output={
                'vcf.gz': '{root}.vcf.gz',
                'viz.gz': '{root}.viz.gz',
                'log.txt': '{root}.log.txt',
            }
        )
        hipstr_job.cloudfuse(f'cpg-{dataset}-main', '/cramfuse')

        # Read in HipSTR variant catalog
        variant_catalog_input = f"{variant_catalog}_shard_{i+1}.bed"
        hipstr_regions = b.read_input(variant_catalog_input)

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
        hipstr_output_path_vcf = output_path(
            f'{output_file_name}_shard{i+1}.vcf.gz', 'analysis'
        )
        b.write_output(hipstr_job.hipstr_output['vcf.gz'], hipstr_output_path_vcf)

        hipstr_output_path_viz = output_path(
            f'{output_file_name}_shard{i+1}.viz.gz', 'analysis'
        )
        b.write_output(hipstr_job.hipstr_output['viz.gz'], hipstr_output_path_viz)

        hipstr_output_path_log = output_path(
            f'{output_file_name}_shard{i+1}.log.txt', 'analysis'
        )
        b.write_output(hipstr_job.hipstr_output['log.txt'], hipstr_output_path_log)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
