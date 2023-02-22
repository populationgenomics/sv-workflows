#!/usr/bin/env python3
# pylint: disable=import-error,duplicate-code

"""
This script uses GangSTRv2.5 to call STRs on WGS cram files.
Required input: --variant-catalog (file path to variant catalog), --dataset, and external sample IDs
For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' str_iterative_gangstr_runner.py --variant-catalog=gs://cpg-tob-wgs-test/hoptan-str/Illuminavariant_catalog.json --dataset=tob-wgs TOB1XXXX TOB1XXXX

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import logging

import click
import hailtop.batch as hb

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.apis import AnalysisApi, SampleApi
from sample_metadata.models import AnalysisStatus

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path, reference_path

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
GANGSTR_IMAGE = config['images']['gangstr']


# inputs:
# variant catalog
@click.option('--variant-catalog', help='Full path to Illumina Variants catalog')
# input dataset
@click.option('--dataset', help='dataset to operate on, eg: tob-wgs')
# input sample ID
@click.argument('external-wgs-ids', nargs=-1)
@click.command()
def main(
    variant_catalog, dataset, external_wgs_ids: list[str]
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=config['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    external_id_to_cpg_id: dict[str, str] = SampleApi().get_sample_id_map_by_external(
        dataset, list(external_wgs_ids)
    )
    cpg_id_to_external_id = {
        cpg_id: external_wgs_id
        for external_wgs_id, cpg_id in external_id_to_cpg_id.items()
    }
    if dataset in ['tob-wgs', 'hgdp']:
        ref_fasta = (
            'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'
        )
    else:
        ref_fasta = reference_path('broad/ref_fasta')
    if dataset == 'tob-wgs':
        analysis_query_model = AnalysisQueryModel(
            sample_ids=list(external_id_to_cpg_id.values()),
            projects=[dataset],
            type=AnalysisType('cram'),
            status=AnalysisStatus('completed'),
            meta={'sequencing_type': 'genome', 'source': 'nagim'},
        )
    else:
        analysis_query_model = AnalysisQueryModel(
            sample_ids=list(external_id_to_cpg_id.values()),
            projects=[dataset],
            type=AnalysisType('cram'),
            status=AnalysisStatus('completed'),
            meta={'sequencing_type': 'genome'},
        )
    crams_path = AnalysisApi().query_analyses(analysis_query_model)
    cpg_sids_with_crams = set(sid for sids in crams_path for sid in sids['sample_ids'])
    cpg_sids_without_crams = set(cpg_id_to_external_id.keys()) - cpg_sids_with_crams
    if cpg_sids_without_crams:
        external_wgs_sids_without_crams = ', '.join(
            cpg_id_to_external_id[sid] for sid in cpg_sids_without_crams
        )
        logging.warning(
            f'There were some samples without CRAMs: {external_wgs_sids_without_crams}'
        )
    gangstr_regions = b.read_input(variant_catalog)

    # Iterate over each sample to call GangSTR
    for cram_obj in crams_path:

        # Making sure Hail Batch would localize both CRAM and the correponding CRAI index
        crams = b.read_input_group(
            **{'cram': cram_obj['output'], 'cram.crai': cram_obj['output'] + '.crai'}
        )
        cpg_sample_id = cram_obj['sample_ids'][0]

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

        # GangSTR job initialisation
        gangstr_job = b.new_job(name=f'GangSTR:{cpg_sample_id} running')
        gangstr_job.image(GANGSTR_IMAGE)
        gangstr_job.storage('50G')
        gangstr_job.cpu(8)

        gangstr_job.declare_resource_group(
            gangstr_output={
                'vcf': '{root}.vcf',
                'insdata': '{root}.insdata.tab',
                'samplestats': '{root}.samplestats.tab',
            }
        )

        gangstr_job.command(
            f"""
        GangSTR --bam {crams['cram']} --ref {ref.base} --regions {gangstr_regions} --out {gangstr_job.gangstr_output} --bam-samps {cpg_sample_id}
        """
        )
        # GangSTR output writing
        gangstr_output_path = output_path(f'{cpg_sample_id}_gangstr')
        b.write_output(gangstr_job.gangstr_output, gangstr_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
