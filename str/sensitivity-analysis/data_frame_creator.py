#!/usr/bin/env python3
# pylint: disable=import-error
"""
analysis-runner --access-level test --dataset hgdp --description 'EH data frame creator' --output-dir 'str/sensitivity-analysis/eh' data_frame_creator.py --input-dir=gs://cpg-hgdp-test/str/sensitivity-analysis/eh

"""
import os
import logging
import pandas as pd



from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path
import hailtop.batch as hb
from google.cloud import storage

config = get_config()


    #csv = ""            


"""
def concatenate_csv(csv_array):
    combo_csv = ""
    for i in csv_array: 
     #   file = 
        combo_csv= combo_csv+i
    return combo_csv"""
def eh_csv_writer(input_dir):
    file = ""
    bucket_name, *components = input_dir[5:].split('/')
    client = storage.Client()
    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    for file in files: 
        if file.endswith(".vcf"): 
            blob = bucket_name.get_blob(input_dir[6+len(bucket_name):])
            with blob.open("r") as f: 
                array = f.readlines() 
                for line in array: 
                    file+= line
    return file

def main():
# pylint: disable=missing-function-docstring
# Initializing Batch
    backend = hb.ServiceBackend(
            billing_project=get_config()['hail']['billing_project'],
            remote_tmpdir=remote_tmpdir(),
        )
    b = hb.Batch(backend= backend, default_python_image=config['workflow']['driver_image'])
    j = b.new_python_job(name = "EH dataframe writer")
    
    #for vcf_file in vcf_path:
    tester = j.call(eh_csv_writer)

    b.write_output(tester.as_str(), output_path('eh_data_frame.txt'))
    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter