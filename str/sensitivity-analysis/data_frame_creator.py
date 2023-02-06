#!/usr/bin/env python3
# pylint: disable=import-error
"""
analysis-runner --access-level test --dataset hgdp --description 'EH data frame creator' --output-dir 'str/sensitivity-analysis/data_frames' data_frame_creator.py --input-dir-eh=gs://cpg-hgdp-test/str/sensitivity-analysis/eh --input-dir-gangstr=gs://cpg-hgdp-test/str/sensitivity-analysis/gangstr

"""
import click

import pandas as pd
from cloudpathlib import AnyPath, CloudPath

from cyvcf2 import VCFReader
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path

config = get_config()


def eh_csv_writer(input_dir):
    #input_dir = 'gs://cpg-hgdp-test/str/sensitivity-analysis/eh'

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

    directory_files = [file for file in AnyPath(input_dir).iterdir() if file.suffix == '.vcf']
    for file in directory_files:

        # use cloudpathlib to copy this file into a local temp file
        tmp_name = f'tmp_{file.name}'
        local_vcf = file.copy(tmp_name)

        # open a vcf reader instance
        for variant in VCFReader(str(local_vcf)):
            chr = variant.CHROM
            start = variant.START
            filter = variant.FILTER or 'PASS'
            # ...

            # getting things from the FORMAT data
            # 'variant_characteristics' below
            # first (and only) index, as this is a single sample VCF
            e_GT = variant.filter('GT'[0])

            # getting things from the INFO string (in dict form)
            csq = variant.INFO.get('CSQ')

            # ...
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
                        )
        )
    return '\n'.join(csv)

def gangstr_csv_writer(input_dir):
   # input_dir = 'gs://cpg-hgdp-test/str/sensitivity-analysis/gangstr' #can't seem to code this as an argument
    bucket_name, *components = input_dir[5:].split('/')
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    csv = (
            '\t'.join(
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
                    start = str(int(attributes[1])-1) #convert back to 0-based to match with EH 
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
                    g_locus_characteristics = attributes[9].split(":")
                    g_GT = g_locus_characteristics[0]
                    g_DP = g_locus_characteristics[1]
                    g_Q = g_locus_characteristics[2]
                    g_REPCN = g_locus_characteristics[3]
                    g_REPCI = g_locus_characteristics[4]
                    g_RC = g_locus_characteristics[5]
                    g_ENCLREADS = g_locus_characteristics[6]
                    g_FLNKREADS = g_locus_characteristics[7]
                    g_ML = g_locus_characteristics[8]
                    g_INS = g_locus_characteristics[9]
                    g_STDERR = g_locus_characteristics[10]
                    g_QEXP = g_locus_characteristics[11]
                    g_allele_1 = g_REPCN.split(",")[0]
                    g_allele_2 = g_REPCN.split(",")[1]
                    csv = csv+(
                        '\t'.join(
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

def merge_csv(eh_csv_obj, gangstr_csv_obj):
    eh_csv = pd.read_csv(eh_csv_obj)
    gangstr_csv = pd.read_csv(gangstr_csv_obj, sep ="\t")
    merged = pd.merge(eh_csv,gangstr_csv, on = ['sample_id', 'chr', 'start'], how = 'outer')
    return merged.to_csv()

@click.command()
@click.option('--input-dir-eh')
@click.option('--input-dir-gangstr')
def main(input_dir_eh, input_dir_gangstr):
# pylint: disable=missing-function-docstring
# Initializing Batch

    eh_csv = eh_csv_writer(AnyPath(input_dir_eh))

    # --- not edited below
    gangstr_csv = gangstr_csv_writer(input_dir_gangstr)
    merged_csv = merge_csv(eh_csv, gangstr_csv)

    b.write_output(eh_csv.as_str(), output_path('eh.csv'))
    b.write_output(gangstr_csv.as_str(), output_path('gangstr.tsv'))
    b.write_output(merged_csv.as_str(), output_path('merged_dataframe.csv'))

    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter