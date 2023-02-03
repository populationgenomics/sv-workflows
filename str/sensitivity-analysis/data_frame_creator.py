import os
import logging

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path, reference_path
from cpg_workflows.batch import get_batch
from google.cloud import storage
config = get_config()

@click.option('--input-dir', help='Full path to directory containing crams')
@click.command()
def main(
    input_dir):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    vcf_path = []
    bucket_name, *components = input_dir[5:].split('/')

    client = storage.Client()

    blobs = client.list_blobs(bucket_name, prefix = '/'.join(components))
    files: Set[str] = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    for file in files: 
        if file.endswith(".vcf"): 
                vcf_path.append(file)

for cram_obj in crams_path:
file= open("str_sensitivity-analysis_eh_CPG19869_eh.vcf")

with open("eh_data_frame.csv", "w", encoding = "utf-8")as handle:
    handle.write(
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
    for line in file: 
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
        handle.write(
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





