#!/usr/bin/env python3
# pylint: disable=import-error

import os
import logging

import click
import hailtop.batch as hb
from cloudpathlib import AnyPath

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.model.analysis_query_model import AnalysisQueryModel
from sample_metadata.apis import AnalysisApi, SampleApi
from sample_metadata.models import AnalysisStatus

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path, reference_path



config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
EH_IMAGE = config['images']['expansionhunter']


# inputs:
@click.command()
def main(
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))


    sample_id_file = "gs://cpg-tob-wgs-test/hoptan-str/karyotype_sex_mapping.csv"
    with AnyPath(sample_id_file).open() as f:
        for line in f:
            index, s, sex, external_id = line.split(",") 
            print(sex)
if __name__ == '__main__':
    main()
