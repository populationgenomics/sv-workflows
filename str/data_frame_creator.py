#!/usr/bin/env python3
# pylint: disable=import-error, too-many-locals, broad-exception-raised
"""
This script merges all the VCFs from one STR caller into a .CSV (ExpansionHunter) or .TSV (GangSTR) format that can be read into R. 
analysis-runner --access-level test --dataset tob-wgs --description 'data frame creator' --output-dir 'hoptan-str/tob_test_crams/data_frames' data_frame_creator.py  --input-dir-eh=gs://cpg-tob-wgs-test/hoptan-str/tob_test_crams/output_calls/eh_0_based --input-dir-gangstr=gs://cpg-tob-wgs-test/hoptan-str/tob_test_crams/output_calls/gangstr_0_based --output-name-eh=eh.csv --output-name-gangstr=gangstr.tsv

"""
import click
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path
import hailtop.batch as hb
from google.cloud import storage

config = get_config()


def eh_csv_writer(input_dir):
    """Creates a CSV file containing dataframe of merged EH VCFs"""
    bucket_name, *components = input_dir[5:].split('/')
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = client.list_blobs(bucket_name, prefix='/'.join(components))
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
                'e_gt',
                'e_so',
                'e_allele_1',
                'e_allele_2',
                'e_repci',
                'e_adsp',
                'e_adfl',
                'e_adir',
                'e_lc',
            ]
        )
        + '\n'
    )

    for file in files:
        if file.endswith('.vcf'):
            blob = bucket.blob(file[6 + len(bucket_name) :])
            with blob.open('r') as f:
                array = f.readlines()
                for line in array:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#CHROM'):
                        header = line.split()
                        sample_id = header[9]
                        continue
                    attributes = line.split()
                    chr = attributes[0]
                    start = attributes[1]
                    e_qual = attributes[6]
                    locus_characteristics = attributes[7].split(';')
                    end = locus_characteristics[0][4:]
                    repeat_units_in_ref = locus_characteristics[1][4:]
                    ref_sequence_length = locus_characteristics[2][3:]
                    motif = locus_characteristics[3][3:]
                    variant_characteristics = attributes[9].split(':')
                    e_gt = variant_characteristics[0]
                    e_so = variant_characteristics[1]
                    e_allele_1 = variant_characteristics[2].split('/')[0]
                    e_allele_2 = variant_characteristics[2].split('/')[1]
                    e_repci = variant_characteristics[3]
                    e_adsp = variant_characteristics[4]
                    e_adfl = variant_characteristics[5]
                    e_adir = variant_characteristics[6]
                    e_lc = variant_characteristics[7]
                    csv = csv + (
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
                                e_gt,
                                e_so,
                                e_allele_1,
                                e_allele_2,
                                e_repci,
                                e_adsp,
                                e_adfl,
                                e_adir,
                                e_lc,
                            ]
                        )
                        + '\n'
                    )
    return csv


def gangstr_csv_writer(input_dir):
    """Creates a TSV file containing dataframe of merged GangSTR VCFs"""
    bucket_name, *components = input_dir[5:].split('/')
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = client.list_blobs(bucket_name, prefix='/'.join(components))
    files = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    csv = (
        '\t'.join(
            [
                'sample_id',
                'chr',
                'start',
                'g_gt',
                'g_dp',
                'g_q',
                'g_repci',
                'g_rc',
                'g_enclreads',
                'g_flnkreads',
                'g_ml',
                'g_ins',
                'g_stderr',
                'g_allele_1',
                'g_allele_2',
                'g_qexp',
            ]
        )
        + '\n'
    )
    for file in files:
        if file.endswith('.vcf'):
            blob = bucket.blob(file[6 + len(bucket_name) :])
            with blob.open('r') as f:
                array = f.readlines()
                for line in array:
                    if line.startswith('##'):
                        continue
                    if line.startswith('#CHROM'):
                        header = line.split()
                        sample_id = header[9]
                        continue
                    attributes = line.split()
                    chr = attributes[0]
                    start = attributes[1]
                    if attributes[9] == '.':  # ie variant is not called
                        g_gt = attributes[9]
                        g_dp = attributes[9]
                        g_q = attributes[9]
                        g_repcn = attributes[9]
                        g_allele_1 = attributes[9]
                        g_allele_2 = attributes[9]
                        g_repci = attributes[9]
                        g_rc = attributes[9]
                        g_enclreads = attributes[9]
                        g_flnkreads = attributes[9]
                        g_ml = attributes[9]
                        g_ins = attributes[9]
                        g_stderr = attributes[9]
                        g_qexp = attributes[9]
                        continue
                    g_locus_characteristics = attributes[9].split(':')
                    g_gt = g_locus_characteristics[0]
                    g_dp = g_locus_characteristics[1]
                    g_q = g_locus_characteristics[2]
                    g_repcn = g_locus_characteristics[3]
                    g_repci = g_locus_characteristics[4]
                    g_rc = g_locus_characteristics[5]
                    g_enclreads = g_locus_characteristics[6]
                    g_flnkreads = g_locus_characteristics[7]
                    g_ml = g_locus_characteristics[8]
                    g_ins = g_locus_characteristics[9]
                    g_stderr = g_locus_characteristics[10]
                    g_qexp = g_locus_characteristics[11]
                    g_allele_1 = g_repcn.split(',')[0]
                    g_allele_2 = g_repcn.split(',')[1]
                    csv = csv + (
                        '\t'.join(
                            [
                                sample_id,
                                chr,
                                start,
                                g_gt,
                                g_dp,
                                g_q,
                                g_repci,
                                g_rc,
                                g_enclreads,
                                g_flnkreads,
                                g_ml,
                                g_ins,
                                g_stderr,
                                g_allele_1,
                                g_allele_2,
                                g_qexp,
                            ]
                        )
                        + '\n'
                    )
    return csv


@click.command()
@click.option('--input-dir-eh', help='Input directory for ExpansionHunter VCFs')
@click.option('--input-dir-gangstr', help='Input directory for GangSTR VCFs')
@click.option('--output-name-eh', help='Output file name for ExpansionHunter eg eh.csv')
@click.option(
    '--output-name-gangstr', help='Output file name for GangSTR eg gangstr.tsv'
)
def main(input_dir_eh, input_dir_gangstr, output_name_eh, output_name_gangstr):
    # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(
        backend=backend, default_python_image=config['workflow']['driver_image']
    )
    j = b.new_python_job(name='EH dataframe writer')
    g = b.new_python_job(name='GangSTR dataframe writer')

    eh_csv = j.call(eh_csv_writer, input_dir_eh)
    gangstr_csv = g.call(gangstr_csv_writer, input_dir_gangstr)

    b.write_output(eh_csv.as_str(), output_path({output_name_eh}))
    b.write_output(gangstr_csv.as_str(), output_path({output_name_gangstr}))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
