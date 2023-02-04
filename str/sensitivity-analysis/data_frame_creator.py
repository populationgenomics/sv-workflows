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

def main():
# pylint: disable=missing-function-docstring
# Initializing Batch
    backend = hb.ServiceBackend(
            billing_project=get_config()['hail']['billing_project'],
            remote_tmpdir=remote_tmpdir(),
        )
    b = hb.Batch(backend= backend, default_python_image=config['workflow']['driver_image'])
    j = b.new_python_job(name = "EH dataframe writer")
    input_dir = "gs://cpg-hgdp-test/str/sensitivity-analysis/eh"
    vcf_path = []
    bucket_name, *components = input_dir[5:].split('/')

    client = storage.Client()

    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    for file in files: 
        if file.endswith(".vcf"): 
                vcf_path.append(file)
    #for vcf_file in vcf_path:
    def eh_csv_writer():
        file= ("gs://cpg-hgdp-test/str/sensitivity-analysis/eh/CPG19869_eh.vcf")
        file = pd.read_csv(file, sep='\t', skiprows = 100)
        for line in file: 
            print(line)
    tester = j.call(eh_csv_writer)

    b.write_output(tester, output_path('eh_data_frame.csv'))
    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter