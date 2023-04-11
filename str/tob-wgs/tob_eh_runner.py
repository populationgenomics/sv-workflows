#!/usr/bin/env python3
# pylint: disable=import-error

"""
This script uses ExpansionHunter v5 to call STRs on the entire cohort of TOB-WGS (1057 samples)
Required input: --variant-catalog (file path to variant catalog), --sample-id-file
For example:
analysis-runner --access-level full --dataset tob-wgs --description 'TOB Run' --output-dir 'tester' tob_eh_runner.py --variant-catalog=gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json --sample-id-file=gs://cpg-tob-wgs-test/hoptan-str/karyotype_sex_mapping.csv

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os

import click
import hailtop.batch as hb

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.apis import AnalysisApi
from sample_metadata.models import AnalysisStatus
from cloudpathlib import AnyPath


from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
EH_IMAGE = config['images']['expansionhunter']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
# sample id and sex mapping file
@click.option(
    '--sample-id-file', help='Full path to mapping of CPG id, TOB id, and sex'
)
@click.command()
def main(variant_catalog, sample_id_file):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    ref_fasta = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
    eh_regions = b.read_input(variant_catalog)

    with AnyPath(sample_id_file).open() as f:

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
                # 'X' and 'ambiguous' karyotypic sex will be marked as female (ExpansionHunter defaults to female if no sex_parameter is provided)

            analysis_query_model = AnalysisQueryModel(
                sample_ids=[cpg_id],
                projects=['tob-wgs'],
                type=AnalysisType('cram'),
                status=AnalysisStatus('completed'),
                meta={'sequencing_type': 'genome', 'source': 'nagim'},
            )

            cram_path = AnalysisApi().query_analyses(analysis_query_model)

            # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
            crams = b.read_input_group(
                **{
                    'cram': cram_path[0]['output'],
                    'cram.crai': cram_path[0]['output'] + '.crai',
                }
            )

            # Working with CRAM files requires the reference fasta
            ref = b.read_input_group(
                **dict(
                    base=ref_fasta,
                    fai=ref_fasta + '.fai',
                    dict=ref_fasta.replace('.fasta', '.dict')
                )
            )

            # ExpansionHunter job initialisation
            eh_job = b.new_job(name=f'ExpansionHunter:{cpg_id}')
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
            --output-prefix {eh_job.eh_output} \\
            --sex {sex_param}
            """
            )
            # ExpansionHunter output writing
            eh_output_path = output_path(f'{cpg_id}_eh')
            b.write_output(eh_job.eh_output, eh_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
