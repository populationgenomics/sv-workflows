#!/usr/bin/env python3


"""
This script uses HipSTR to call STRs on WGS cram files, using the joint calling option, and a sharded catalog.
For example:
analysis-runner --access-level test --dataset tob-wgs --description 'hipstr run' --output-dir 'str/sensitivity-analysis/hipstr' str_iterative_hipstr_runner.py  --output-file-name=hipster_90_genomes  --variant-catalog=gs://.... --dataset=hgdp --max-str-len=100 --sample-id-file=gs://...

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import logging

import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path
from metamist.graphql import gql, query

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
        """,
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
        logging.warning(f'There were some samples without CRAMs: {cpg_sids_without_crams}')

    # Create string containing paths based on /cramfuse
    cramfuse_path = []
    for cram_obj in crams_by_id.values():
        suffix = cram_obj['output'].removeprefix('gs://').split('/', maxsplit=1)[1]
        cramfuse_path.append(f'/cramfuse/{suffix}')
    cramfuse_path = ','.join(cramfuse_path)  # string format for input into hipstr
    return cramfuse_path


# inputs:
@click.option('--sample-id-file', help='Path to CSV of CPG IDs')
@click.option('--job-storage', help='Storage of the Hail batch job eg 30G', default='30G')
@click.option('--job-memory', help='Memory of the Hail batch job', default='highmem')
@click.option(
    '--variant-catalog',
    help='Full path to HipSTR Variants sharded catalog directory or file path to unsharded catalog',
)
@click.option('--dataset', help='dataset eg tob-wgs')
@click.option('--output-file-name', help='Output file name without file extension')
@click.option(
    '--max-str-len',
    help="Only genotype STRs in the provided BED file with length < MAX_BP (Default = 100)",
    default='100',
)
@click.command()
def main(
    sample_id_file,
    job_storage,
    job_memory,
    variant_catalog,
    dataset,
    output_file_name,
    max_str_len,
):  # pylint: disable=missing-function-docstring
    b = get_batch()
    ref_fasta = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
    # Read in reference
    ref = b.read_input_group(
        **dict(
            base=ref_fasta,
            fai=ref_fasta + '.fai',
            dict=ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict',
        ),
    )
    internal_cpg_ids = []
    # open sample id file and add to list
    with to_path(sample_id_file).open() as f:
        for line in f:
            split_line = line.split(',')
            cpg_id = split_line[0]
            if cpg_id == 's':  # header line
                continue
            internal_cpg_ids.append(cpg_id)

    catalog_files = list(to_path(variant_catalog).glob('*.bed'))

    cramfuse_path = get_cloudfuse_paths(dataset, internal_cpg_ids)

    for index, subcatalog in enumerate(catalog_files):
        # Create HipSTR job
        hipstr_job = b.new_job(name=f'HipSTR shard {index+1} running')
        hipstr_job.image(HIPSTR_IMAGE)
        hipstr_job.storage(job_storage)
        hipstr_job.cpu(4)
        hipstr_job.memory(job_memory)

        hipstr_job.declare_resource_group(
            hipstr_output={
                'vcf.gz': '{root}.vcf.gz',
                'viz.gz': '{root}.viz.gz',
                'log.txt': '{root}.log.txt',
            },
        )
        hipstr_job.cloudfuse(f'cpg-{dataset}-main', '/cramfuse')

        # Read in HipSTR variant catalog

        hipstr_regions = b.read_input(subcatalog)

        hipstr_job.command(
            f"""
        HipSTR --bams {cramfuse_path} \\
            --fasta {ref.base} \\
            --regions {hipstr_regions} \\
            --str-vcf {hipstr_job.hipstr_output['vcf.gz']} \\
            --viz-out {hipstr_job.hipstr_output['viz.gz']} \\
            --log {hipstr_job.hipstr_output['log.txt']} \\
            --max-str-len {max_str_len} \\
            --output-filters
        """,
        )
        # HipSTR output writing
        hipstr_output_path_name = output_path(f'{output_file_name}_shard{index+1}', 'analysis')
        b.write_output(hipstr_job.hipstr_output, hipstr_output_path_name)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
