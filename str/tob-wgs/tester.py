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

def main():
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    with AnyPath('gs://cpg-tob-wgs-test/hoptan-str/karyotype_sex_mapping.csv').open() as f:


        # Iterate over each sample to call Expansion Hunter
        for line in f:
            cpg_id, sex, external_id = line.split(",")
            if cpg_id == "s":  # header line
                continue
            print(sex)
            if sex == "XY":
                sex_param = "male"
            else:
                sex_param = "female"
                # "X" and "ambiguous" karyotypic sex will be marked as female (ExpansionHunter defaults to female if no sex_parameter is provided)

            analysis_query_model = AnalysisQueryModel(
                sample_ids=list([cpg_id]),
                projects=['tob-wgs'],
                type=AnalysisType('cram'),
                status=AnalysisStatus('completed'),
                meta={'sequencing_type': 'genome', 'source': 'nagim'},
            )

            cram_path = AnalysisApi().query_analyses(analysis_query_model)
            potato = cram_path[0]['output']
            print(f'{potato} {sex_param}')
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter

