#!/usr/bin/env python3
# pylint: disable=import-error

import hailtop.batch as hb

from cpg_utils.config import get_config

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
EH_IMAGE = config['images']['expansionhunter']


# inputs:
@click.command()
def main(
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))


    sample_id_file = "gs://cpg-tob-wgs-test/hoptan-str/karyotype_sex_mapping.csv"
    pointer = b.read_input(sample_id_file)
    with open (pointer) as f:
        for line in f: 
            print(line)
