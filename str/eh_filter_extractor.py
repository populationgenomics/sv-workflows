#!/usr/bin/env python3
# pylint: disable=import-error, too-many-locals

"""
This script extracts the EH binary filter (PASS/LowDepth) for all EH VCFs found within a directory into a .TSV format to be used as an annotation file in Hail query. 

analysis-runner --access-level standard --dataset tob-wgs --description \
    'EH Filter Extractor' --output-dir 'str/expansionhunter/v3-eh-filter-extractor' \
    eh_filter_extractor.py \
    --input-dir-eh=gs://cpg-tob-wgs-main/str/expansionhunter/v3 \
    --output-name-eh=eh.tsv 
"""


import click
from cloudpathlib import GSPath
from cyvcf2 import VCFReader
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path


config = get_config()


def eh_filter_extractor(input_dir):
    """Creates a TSV file containing locus, sample_id, and binary filter status of EH VCFs"""
    files = to_path(input_dir).glob('*.vcf')
    tsv = (
        '\t'.join(
            [
                'sample_id',
                'locus',
                'e_qual',
            ]
        )
        + '\n'
    )

    for file in files:
        if isinstance(file, GSPath):
            file.copy(file.name)
            file = file.name
        reader = VCFReader(file)
        sample_id = str(reader.samples[0])
        for variant in reader:
            chr = str(variant.CHROM)
            start = str(variant.POS)
            locus = f'{chr}:{start}'
            e_qual = str(variant.FILTER)
            tsv = tsv + ('\t'.join([sample_id, locus, e_qual]) + '\n')
    return tsv


@click.command()
@click.option('--input-dir-eh', help='Input directory for ExpansionHunter VCFs')
@click.option(
    '--output-name-eh',
    help='Output file name for filter extractor output file eg eh.tsv',
)
def main(
    input_dir_eh,
    output_name_eh,
):
    """
    some docstring
    """
    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(
        backend=backend, default_python_image=config['workflow']['driver_image']
    )
    j = b.new_python_job(name='EH filter status extractor')

    eh_tsv = j.call(eh_filter_extractor, input_dir_eh)

    b.write_output(eh_tsv.as_str(), output_path(output_name_eh, 'analysis'))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
