#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses ExpansionHunterv5 to call STRs on WGS cram files.
Required input: --variant-catalog (file path to variant catalog), --dataset, and internal sample IDs
For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' str_iterative_eh_runner.py --variant-catalog=gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json --dataset=tob-wgs CPG26 CPG18

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import logging

import click
import hailtop.batch as hb

from metamist.graphql import gql, query

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path, reference_path

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
EH_IMAGE = config['images']['expansionhunter']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
# input dataset
@click.option('--dataset', help='dataset eg tob-wgs')
# input sample ID
@click.argument('input_cpg_sids', nargs=-1)
@click.command()
def main(
    variant_catalog, dataset, input_cpg_sids: list[str]
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    if dataset in ['tob-wgs', 'hgdp']:
        ref_fasta = (
            'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
        )
    else:
        ref_fasta = reference_path('broad/ref_fasta')

    cram_retrieval_query = gql(
        """
        query MyQuery($dataset: String!,$input_cpg_sids: [String!]!) {
    project(name: $dataset) {
        sequencingGroups(id: {in_: $input_cpg_sids}) {
            id
            sample {
                externalId
            }
            analyses(type: {eq: 'cram'}, active: {eq: true}) {
                output
                timestampCompleted
            }
        }
    }
        """
    )
    response = query(
        cram_retrieval_query,
        variables={'dataset': dataset, 'input_cpg_sids': input_cpg_sids},
    )
    valid_cpg_ids = []
    crams_path = []
    for i in response['project']['sequencingGroups']:
        valid_cpg_ids.append(i['id'])
        if (
            len(i['analyses']) > 1
        ):  # if there are multiple crams stored per CPG sample ID
            for cram in i['analyses']:
                if (
                    'archive' not in cram['output']
                ):  # select the cram that is not archived
                    cram['id'] = i['id']  # ensures cpg id is stored in the dictionary
                    crams_path.append(cram)
                    break

    if len(valid_cpg_ids) != len(input_cpg_sids):
        cpg_sids_without_crams = set(input_cpg_sids) - set(valid_cpg_ids)
        logging.warning(
            f'There were some samples without CRAMs: {cpg_sids_without_crams}'
        )

    eh_regions = b.read_input(variant_catalog)

    # Iterate over each sample to call Expansion Hunter
    for cram_obj in crams_path:
        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj['output'], 'cram.crai': cram_obj['output'] + '.crai'}
        )
        cpg_sample_id = cram_obj['id']

        # Working with CRAM files requires the reference fasta
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

        # ExpansionHunter job initialisation
        eh_job = b.new_job(name=f'ExpansionHunter:{cpg_sample_id} running')
        eh_job.image(EH_IMAGE)
        eh_job.storage('50G')
        eh_job.cpu(8)

        eh_job.declare_resource_group(
            eh_output={
                'vcf': '{root}.vcf',
                'json': '{root}.json',
                'realigned.bam': '{root}_realigned.bam',
            }
        )

        eh_job.command(
            f"""
        ExpansionHunter  \\
        --reads {crams['cram']} \\
        --reference {ref.base} --variant-catalog {eh_regions}\\
         --threads 16 --analysis-mode streaming \\
         --output-prefix {eh_job.eh_output}
        """
        )
        # ExpansionHunter output writing
        eh_output_path = output_path(f'{cpg_sample_id}_eh', 'analysis')
        b.write_output(eh_job.eh_output, eh_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
