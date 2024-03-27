#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member

"""
Script to submit a dataproc job

"""
from analysis_runner import dataproc
from cpg_utils.hail_batch import get_batch

def main():
    script = (
        f'dataproc_helper/str_pca_dataproc_hail_script.py '
        f'--file-path=gs://cpg-bioheart-test/str/associatr/mt_filtered/v1/str.mt'
    )
    j = dataproc.hail_dataproc_job(
                get_batch(),
                script,
                max_age='48h',
                packages=[
                    'cpg_workflows',
                    'google',
                    'gcloud',
                ],
                num_workers=2,
                num_secondary_workers=0,
                job_name='STR-PCA',
            )
    j._preemptible = False
    j.attributes = (j.attributes or {}) | {'tool': 'hailctl dataproc'}
    get_batch.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter