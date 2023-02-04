#!/usr/bin/env python3
# pylint: disable=import-error
"""
analysis-runner --access-level test --dataset hgdp --description 'EH data frame creator' --output-dir 'str/sensitivity-analysis/data_frames' data_frame_creator.py --input-dir=gs://cpg-hgdp-test/str/sensitivity-analysis/eh

"""
import os
import logging
import click
import pandas as pd


from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path
import hailtop.batch as hb
from google.cloud import storage

config = get_config()


def eh_csv_writer():
    input_dir = 'gs://cpg-hgdp-test/str/sensitivity-analysis/eh' #can't seem to code this as an argument
    bucket_name, *components = input_dir[5:].split('/')
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    csv = (
            ','.join(
                [   
                    'sample_id',
                    'chr',
                    'start',
                    'e_qual',
                    'end',
                    'repeat_units_in_ref',
                    'ref_sequence_length',
                    'motif',
                    'e_GT',
                    'e_SO',
                    'e_allele_1',
                    'e_alelle_2',
                    'e_REPCI',
                    'e_ADSP',
                    'e_ADFL',
                    'e_ADIR',
                    'e_LC'
                ]
            )
            +'\n'
        )

    for file in files: 
        if file.endswith(".vcf"): 
            blob = bucket.blob(file[6+len(bucket_name):])
            with blob.open("r") as f: 
                array = f.readlines() 
                for line in array:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#CHROM"):
                        header = line.split()
                        sample_id = header[9]
                        continue
                    attributes = line.split()
                    chr = attributes[0]
                    start = attributes[1]
                    e_qual = attributes[6]
                    locus_characteristics = attributes[7].split(";")
                    end = locus_characteristics[0][4:]
                    repeat_units_in_ref = locus_characteristics[1][4:]
                    ref_sequence_length = locus_characteristics[2][3:]
                    motif = locus_characteristics[3][3:]
                    variant_characteristics = attributes[9].split(":")
                    e_GT = variant_characteristics[0]
                    e_SO = variant_characteristics[1]
                    e_allele_1 = variant_characteristics[2].split("/")[0]
                    e_allele_2 = variant_characteristics[2].split("/")[1]
                    e_REPCI = variant_characteristics[3]
                    e_ADSP = variant_characteristics[4]
                    e_ADFL = variant_characteristics[5]
                    e_ADIR = variant_characteristics[6]
                    e_LC = variant_characteristics[7]
                    csv = csv+(
                        ','.join(
                            [
                                sample_id, 
                                chr, 
                                start, 
                                e_qual, 
                                end,
                                repeat_units_in_ref,
                                ref_sequence_length, 
                                motif,
                                e_GT,
                                e_SO, 
                                e_allele_1,
                                e_allele_2,
                                e_REPCI,
                                e_ADSP, 
                                e_ADFL,
                                e_ADIR,
                                e_LC
                            ]
                        )+'\n'
        )
    return csv

def gangstr_csv_writer():
    input_dir = 'gs://cpg-hgdp-test/str/sensitivity-analysis/gangstr' #can't seem to code this as an argument
    bucket_name, *components = input_dir[5:].split('/')
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    csv = (
            ','.join(
                [   
                    'sample_id',
                    'chr',
                    'start',
                    'g_GT',
                    'g_DP',
                    'g_Q',
                    'g_REPCI',
                    'g_RC',
                    'g_ENCLREADS',
                    'g_FLNKREADS',
                    'g_ML',
                    'g_INS',
                    'g_STDERR',
                    'g_allele_1',
                    'g_allele_2',
                    'g_QEXP'
                ]
            )
            +'\n'
        )
    for file in files: 
        if file.endswith(".vcf"): 
            blob = bucket.blob(file[6+len(bucket_name):])
            with blob.open("r") as f: 
                array = f.readlines() 
                for line in array:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#CHROM"):
                        header = line.split()
                        sample_id = header[9]
                        continue
                    attributes = line.split()
                    chr = attributes[0]
                    start = attributes[1]
                    if attributes[9] == ".": #ie variant is not called
                        g_GT = attributes[9]
                        g_DP = attributes[9]
                        g_Q = attributes[9]
                        g_REPCN = attributes[9]
                        g_allele_1 = attributes[9]
                        g_allele_2=attributes[9]
                        g_REPCI = attributes[9]
                        g_RC = attributes[9]
                        g_ENCLREADS = attributes[9]
                        g_FLNKREADS = attributes[9]
                        g_ML = attributes[9]
                        g_INS = attributes[9]
                        g_STDERR = attributes[9]
                        g_QEXP = attributes[9]
                        continue
                    locus_characteristics = attributes[9].split(":")
                    g_GT = locus_characteristics[0]
                    g_DP = locus_characteristics[1]
                    g_Q = locus_characteristics[2]
                    g_REPCN = locus_characteristics[3]
                    g_allele_1 = g_REPCN.split(",")[0]
                    g_allele_2 = g_REPCN.split(",")[1]
                    g_REPCI = locus_characteristics[4]
                    g_RC = locus_characteristics[5]
                    g_ENCLREADS = locus_characteristics[6]
                    g_FLNKREADS = locus_characteristics[7]
                    g_ML = locus_characteristics[8]
                    g_INS = locus_characteristics[9]
                    g_STDERR = locus_characteristics[10]
                    g_QEXP = locus_characteristics[11]
                    csv = csv+(
                        ','.join(
                            [
                                sample_id, 
                                chr, 
                                start, 
                                g_GT, 
                                g_DP,
                                g_Q,
                                g_REPCI, 
                                g_RC,
                                g_ENCLREADS,
                                g_FLNKREADS, 
                                g_ML,
                                g_INS,
                                g_STDERR,
                                g_allele_1, 
                                g_allele_2,
                                g_QEXP
                            ]
                        )+'\n'
        )
    return csv

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

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter