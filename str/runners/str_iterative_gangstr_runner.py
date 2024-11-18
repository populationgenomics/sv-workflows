#!/usr/bin/env python3
# pylint: disable=import-error,duplicate-code

"""
This script uses GangSTRv2.5 to call STRs on WGS cram files.
Required input: --variant-catalog (file path to variant catalog), --dataset, and external sample IDs
For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' str_iterative_gangstr_runner.py --variant-catalog=gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json --dataset=tob-wgs-test CPG308239

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import logging
import os

import click


import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, reference_path

from cpg_utils.hail_batch import get_batch, output_path, reference_path
from metamist.graphql import gql, query

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
GANGSTR_IMAGE = config['images']['gangstr']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
# input dataset
@click.option('--dataset', help='dataset to operate on, eg: tob-wgs')
# input sample ID
@click.argument('cpg-ids', nargs=-1)
@click.command()
def main(variant_catalog, dataset, cpg_ids: list[str]):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    if any(keyword in dataset for keyword in ['hgdp', 'thousand-genomes']):
        ref_fasta = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
    else:
        ref_fasta = reference_path('broad/ref_fasta')

    gangstr_regions = b.read_input(variant_catalog)

    # Working with CRAM files requires the reference fasta
    ref = b.read_input_group(
        **dict(
            base=ref_fasta,
            fai=ref_fasta + '.fai',
            dict=ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict',
        ),
    )

    # Iterate over each sample to call GangSTR
    for cpg_id in cpg_ids:

         # retrieve corresponding cram path
        cram_retrieval_query = gql(
            """
            query MyQuery($dataset: String!,$cpg_id: [String!]!) {
            project(name: $dataset) {
            sequencingGroups(id: {in_: $cpg_id}) {
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
            variables={'dataset': dataset, 'cpg_id': cpg_id},
        )

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{
                'cram': response['project']['sequencingGroups'][0]['analyses'][0]['output'],
                'cram.crai': response['project']['sequencingGroups'][0]['analyses'][0]['output'] + '.crai',
            },
        )

        # GangSTR job initialisation
        gangstr_job = b.new_job(name=f'GangSTR:{cpg_id} running')
        gangstr_job.image(GANGSTR_IMAGE)
        gangstr_job.storage('50G')
        gangstr_job.cpu(8)

        gangstr_job.declare_resource_group(
            gangstr_output={
                'vcf': '{root}.vcf',
                'insdata': '{root}.insdata.tab',
                'samplestats': '{root}.samplestats.tab',
            },
        )

        gangstr_job.command(
            f"""
        GangSTR --bam {crams['cram']} --ref {ref.base} --regions {gangstr_regions} --out {gangstr_job.gangstr_output} --bam-samps {cpg_id}
        """,
        )
        # GangSTR output writing
        gangstr_output_path = output_path(f'{cpg_id}_gangstr')
        b.write_output(gangstr_job.gangstr_output, gangstr_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
