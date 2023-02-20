#!/usr/bin/env python3
# pylint: disable=import-error, too-many-locals, broad-exception-raised
#!export HTSLIB_CONFIGURE_OPTIONS="--enable-gcs"
#!apt-get install libbz2-dev libcurl4-openssl-dev
#!pip3 install pysam -v --force-reinstall --no-binary :all:
"""
This script merges all the VCFs from one STR caller into a .CSV (ExpansionHunter) or .TSV (GangSTR) format that can be read into R. 
analysis-runner --access-level test --dataset tob-wgs --description 'data frame creator' --output-dir 'hoptan-str/tob_test_crams/data_frames' data_frame_creator.py  --input-dir-eh=gs://cpg-tob-wgs-test/hoptan-str/tob_test_crams/output_calls/eh_0_based --input-dir-gangstr=gs://cpg-tob-wgs-test/hoptan-str/tob_test_crams/output_calls/gangstr_0_based --output-name-eh=eh.csv --output-name-gangstr=gangstr.tsv

"""
from cyvcf2 import VCF
import click
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path
import hailtop.batch as hb
from cpg_utils import to_path
from google.cloud import storage
import os

config = get_config()


def eh_csv_writer(input_dir):
    """Creates a CSV file containing dataframe of merged EH VCFs"""
    bucket_name, *components = input_dir[5:].split('/')
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blobs = client.list_blobs(bucket_name, prefix='/'.join(components))
    files = {f'gs://{bucket_name}/{blob.name}' for blob in blobs}
    csv = (
        ",".join(
            [
                "sample_id",
                "chr",
                "start",
                "e_qual",
                "end",
                "repeat_units_in_ref",
                "ref_sequence_length",
                "motif",
                "e_gt",
                "e_so",
                "e_allele_1",
                "e_allele_2",
                "e_repci",
                "e_adsp",
                "e_adfl",
                "e_adir",
                "e_lc",
            ]
        )
        + "\n"
    )

    for file_path in files:
        if file_path.endswith('.vcf'):
            blob = bucket.blob(file_path[6 + len(bucket_name) :])
            with blob.open('r') as file:
                sample_id = str(VCF(file).samples[0])
                for variant in VCF(file):
                    chr = str(variant.CHROM)
                    start = str(variant.POS)
                    e_qual = str(variant.FILTER)
                    end = str(variant.INFO.get("END"))
                    repeat_units_in_ref = str(variant.INFO.get("REF"))
                    ref_sequence_length = str(variant.INFO.get("RL"))
                    motif = str(variant.INFO.get("RU"))
                    e_gt = f"{variant.genotypes[0][0]}/{variant.genotypes[0][1]}"
                    e_so = str(variant.format("SO")[0])
                    e_allele_1 = str(variant.format("REPCN")[0].split("/")[0])
                    e_allele_2 = str(variant.format("REPCN")[0].split("/")[1])
                    e_repci = str(variant.format("REPCI")[0])
                    e_adsp = str(variant.format("ADSP")[0])
                    e_adfl = str(variant.format("ADFL")[0])
                    e_adir = str(variant.format("ADIR")[0])
                    e_lc = str(variant.format("LC")[0][0])
                    csv = csv + (
                        ",".join(
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
                        + "\n"
                    )
    return csv


def gangstr_tsv_writer(input_dir):
    """Creates a TSV file containing dataframe of merged GangSTR VCFs"""
    files = to_path(input_dir).glob("*.vcf")
    csv = (
        "\t".join(
            [
                "sample_id",
                "chr",
                "start",
                "g_gt",
                "g_dp",
                "g_q",
                "g_repci",
                "g_rc",
                "g_enclreads",
                "g_flnkreads",
                "g_ml",
                "g_ins",
                "g_stderr",
                "g_allele_1",
                "g_allele_2",
                "g_qexp",
            ]
        )
        + "\n"
    )
    for file in files:
        sample_id = str(VCF(file).samples[0])
        for variant in VCF(file):
            chr = str(variant.CHROM)
            start = str(variant.POS)
            g_gt = f"{variant.genotypes[0][0]}/{variant.genotypes[0][1]}"
            g_dp = str(variant.format("DP")[0][0])
            g_q = str(variant.format("Q")[0][0])
            g_repcn = (
                f'{variant.format("REPCN")[0][0]}, {variant.format("REPCN")[0][1]}'
            )
            g_repci = str(variant.format("REPCI")[0])
            g_rc = str(variant.format("RC")[0])
            g_enclreads =str( variant.format("ENCLREADS")[0])
            g_flnkreads = str(variant.format("FLNKREADS")[0])
            g_ml = str(variant.format("ML")[0][0])
            g_ins = f'{variant.format("INS")[0][0]},{variant.format("INS")[0][1]}'
            g_stderr = (
                f'{variant.format("STDERR")[0][0]},{variant.format("STDERR")[0][1]}'
            )
            g_qexp = f'{variant.format("QEXP")[0][0]},{variant.format("QEXP")[0][1]},{variant.format("QEXP")[0][2]}'
            g_allele_1 = str(g_repcn.split(",")[0].strip())
            g_allele_2 = str(g_repcn.split(",")[1].strip())
            csv = csv + (
                "\t".join(
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
                + "\n"
            )
    return csv


@click.command()
@click.option("--input-dir-eh", help="Input directory for ExpansionHunter VCFs")
@click.option("--input-dir-gangstr", help="Input directory for GangSTR VCFs")
@click.option("--output-name-eh", help="Output file name for ExpansionHunter eg eh.csv")
@click.option(
    "--output-name-gangstr", help="Output file name for GangSTR eg gangstr.tsv"
)
def main(input_dir_eh, input_dir_gangstr, output_name_eh, output_name_gangstr):
    # pylint: disable=missing-function-docstring
    # Initializing Batch
    os.environ['HTSLIB_CONFIGURE_OPTIONS'] = "--enable-gcs"
    backend = hb.ServiceBackend(
        billing_project=get_config()["hail"]["billing_project"],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(
        backend=backend, default_python_image=config["workflow"]["driver_image"]
    )
    j = b.new_python_job(name="EH dataframe writer")
    g = b.new_python_job(name="GangSTR dataframe writer")

    eh_csv = j.call(eh_csv_writer, input_dir_eh)
    gangstr_tsv = g.call(gangstr_tsv_writer, input_dir_gangstr)

    b.write_output(eh_csv.as_str(), output_path(output_name_eh, "analysis"))
    b.write_output(gangstr_tsv.as_str(), output_path(output_name_gangstr, "analysis"))

    b.run(wait=False)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
