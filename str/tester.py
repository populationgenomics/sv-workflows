#!/usr/bin/env python3

import os
import hailtop.batch as hb
import click

from sample_metadata.models.analysis_type import AnalysisType
from sample_metadata.models.analysis_query_model import AnalysisQueryModel
from sample_metadata.apis import AnalysisApi

from cpg_utils.config import get_config

from cpg_utils.hail_batch import (
    authenticate_cloud_credentials_in_job,
    copy_common_env,
    remote_tmpdir,
)

config = get_config()

DATASET = config['workflow']['dataset']
HAIL_BUCKET = config['workflow']['hail_bucket']
OUTPUT_SUFFIX = config['workflow']['output_prefix']
BILLING_PROJECT = config['hail']['billing_project']
REF_FASTA = os.path.join(config['workflow']['reference_prefix'], 'hg38/v0/Homo_sapiens_assembly38.fasta')
SAMTOOLS_IMAGE = os.path.join(config['workflow']['image_registry_prefix'], 'samtools:v0')
EH_IMAGE = os.path.join(config['workflow']['image_registry_prefix'], 'expansionhunter:5.0.0')

#required input

#variant catalog
@click.option('--ehregions', help='Full path to Illumina Variants catalog')

#input TOB ID 
@click.option('--id', help = 'TOB-WGS ID eg TOB20000')



@click.command()

def main(ehregions, id):  # pylint: disable=missing-function-docstring 
   # Initializing Batch
   backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    analysis_query_model = b.AnalysisQueryModel(
        sample_ids=["CPG01"],
        projects=[DATASET],
        type=AnalysisType("cram"),
        status=AnalysisStatus("completed"),
    )
