#!/usr/bin/env python3
# pylint: disable=import-error

import os
import logging
import click
import pandas as pd


from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path
import hailtop.batch as hb
from google.cloud import storage

config = get_config()

def merge_csv(eh_csv, gangstr_csv):
    eh_csv = pd.read_csv(eh_csv)
    gangstr_csv = pd.read_csv("gs://cpg-hgdp-test/str/sensitivity-analysis/eh/gangstr.csv")
    return merged



@click.command()
@click.option('--input-dir')
def main(input_dir):
# pylint: disable=missing-function-docstring
# Initializing Batch
    backend = hb.ServiceBackend(
            billing_project=get_config()['hail']['billing_project'],
            remote_tmpdir=remote_tmpdir(),
        )
    b = hb.Batch(backend= backend, default_python_image=config['workflow']['driver_image'])
    j = b.new_python_job(name = "EH dataframe writer")
    g = b.new_python_job(name = "GangSTR dataframe writer")
    combo = b.new_python_job(name ="Combined dataframe writer")
    
    eh_csv = j.call(eh_csv_writer)
    gangstr_csv = g.call(gangstr_csv_writer)
    merged_csv = combo.call(merge_csv(eh_csv.as_str(), gangstr_csv.as_str()))

    b.write_output(eh_csv.as_str(), output_path('eh.csv'))
    b.write_output(gangstr_csv.as_str(), output_path('gangstr.csv'))
    b.write_output(merged_csv.as_str(), output_path('merged.csv'))

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter