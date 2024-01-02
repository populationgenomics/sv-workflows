#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses ExpansionHunterv5 to call STRs on WGS cram files.
Required input: --variant-catalog (file path to variant catalog, can be sharded or unsharded), --dataset, and sample mapping file [TSV] (CPG sample id (first column) and sex (second column))
EH will run on every sample listed inthe sample mapping file.

For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' str_iterative_eh_runner.py --variant-catalog=gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json --dataset=tob-wgs --sample-id-file=gs://...

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os

import click
import hailtop.batch as hb

from metamist.graphql import gql, query

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path, reference_path

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
EH_IMAGE = config['images']['expansionhunter_bw2']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
# input dataset
@click.option('--dataset', help='dataset eg tob-wgs')
@click.option(
    '--max-parallel-jobs',
    type=int,
    default=50,
    help=('To avoid exceeding Google Cloud quotas, set this concurrency as a limit.'),
)
# sample id and sex mapping file
@click.option(
    '--sample-id-file', help='Full path to mapping of CPG id, TOB id, and sex'
)
@click.command()
def main(
    variant_catalog, dataset, max_parallel_jobs, sample_id_file
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    ref_fasta = reference_path('broad/ref_fasta')

    # Working with CRAM files requires the reference fasta
    ref = b.read_input_group(
        **dict(
            base=ref_fasta,
            fai=ref_fasta + '.fai',
            dict=ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )
    catalog_files = list(to_path(variant_catalog).glob('*.bed'))
    jobs = []

    with to_path(sample_id_file).open() as f:
        # Iterate over each sample to call Expansion Hunter
        for line in f:
            cpg_id = line.split(',')[0]
            sex = line.split(',')[1]
            if cpg_id == 's':  # header line
                continue
            if sex == 'XY':
                sex_param = 'male'
            else:
                sex_param = 'female'
                # 'X' and 'ambiguous' karyotypic sex will be marked
                # as female (ExpansionHunter defaults to female if
                # no sex_parameter is provided)

            cram_retrieval_query = gql(
                """
                query MyQuery($dataset: String!,$cpg_id: [String!]) {
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
                """
            )
            response = query(
                cram_retrieval_query,
                variables={'dataset': dataset, 'cpg_id': cpg_id},
            )

            # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
            crams = b.read_input_group(
                **{
                    'cram': response['project']['sequencingGroups'][0]['analyses'][0][
                        'output'
                    ],
                    'cram.crai': response['project']['sequencingGroups'][0]['analyses'][
                        0
                    ]['output']
                    + '.crai',
                }
            )

            for index, subcatalog in enumerate(catalog_files):
                # ExpansionHunter job initialisation
                eh_job = b.new_job(
                    name=f'ExpansionHunter:{cpg_id} running  shard {index+1}/{len(catalog_files)}'
                )
                eh_job.image(EH_IMAGE)
                # limit parallelisation
                if len(jobs) >= max_parallel_jobs:
                    eh_job.depends_on(jobs[-max_parallel_jobs])
                jobs.append(eh_job)
                eh_job.storage('70G')
                eh_job.memory('120G')
                eh_job.cpu(16)
                eh_regions = b.read_input(subcatalog)

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
                --output-prefix {eh_job.eh_output} \\
                --sex {sex_param}
                """
                )
                # ExpansionHunter output writing
                eh_output_path = output_path(f'{cpg_id}_eh_shard{index+1}', 'analysis')
                b.write_output(eh_job.eh_output, eh_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
