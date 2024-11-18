#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda
"""
This script exports the rows of an MT as a TSV file

 analysis-runner --dataset "bioheart" \
    --description "rows exporter" \
    --access-level "test" \
    --output-dir "str/polymorphic_run/mt/bioheart_tob/v1_n2412" \
    rows_exporter.py --mt-path=gs://cpg-bioheart-test/str/polymorphic_run/mt/bioheart_tob/v1_n2412/str_annotated.mt


"""

import click


from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


def qc_filter(mt_path):
    """
    Applies QC filters to the input MT
    """
    import hail as hl
    from cpg_utils.hail_batch import init_batch, output_path
    init_batch(worker_memory='highmem')

    # read in mt
    mt = hl.read_matrix_table(mt_path)


    # write out mt
    mt_path = output_path('str_annotated_rows.tsv.bgz')
    mt.rows().export(mt_path)



@click.option(
    '--mt-path',
    help='GCS file path to mt (output of qc_annotator.py)',
    type=str,
)
@click.option('--hail-storage', help='Hail storage', type=str, default='0G')
@click.option('--hail-cpu', help='Hail CPU', type=int, default=1)
@click.option('--hail-memory', help='Hail memory', type=str, default='standard')
@click.command()
def main(
    mt_path,
    hail_storage,
    hail_cpu,
    hail_memory,
):
    """
    Runner to apply QC filters to input MT, and bgzip and tabix.
    """

    b = get_batch('Apply QC filters to MT and export as MT')

    hail_job = b.new_python_job(name='QC filters')
    hail_job.image(config['workflow']['driver_image'])
    hail_job.storage(hail_storage)
    hail_job.cpu(hail_cpu)
    hail_job.memory(hail_memory)
    hail_job.call(qc_filter, mt_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter