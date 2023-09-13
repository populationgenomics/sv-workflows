#!/usr/bin/env python3

"""
This script extracts sample/REPID combinations with EH binary filter set to "LowDepth" for all EH VCFs found within a directory into a .JSON format to be used as an annotation file in Hail query. 

analysis-runner --access-level standard --dataset tob-wgs --description \
    'EH Filter Extractor' --output-dir 'str/expansionhunter/v3-eh-filter-extractor' \
    eh_filter_extractor.py \
    --input-dir-eh=gs://cpg-tob-wgs-main/str/expansionhunter/v3 \
    --output-name-eh=eh.json 
"""

import json
import logging
from collections import defaultdict
import click
from cloudpathlib import GSPath
from cyvcf2 import VCFReader
import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path


config = get_config()


def eh_filter_extractor(input_dir):
    """Creates a JSON file keyed by REPID with value being a list of sample_ids with binary filter status == "LowDepth" in EH VCFs"""

    files = to_path(input_dir).glob('*.vcf')
    low_depth_dict = defaultdict(list)
    local_file = 'local.vcf'
    for file in files:
        if not isinstance(file, GSPath):
            logging.warning(f'There is a file path that is not a GSPath: {file}')
            continue

        file.copy(local_file)
        reader = VCFReader(local_file)
        sample_id = str(reader.samples[0])
        for variant in reader:
            repid = variant.INFO.get('REPID')
            e_qual = str(variant.FILTER)
            if e_qual == 'LowDepth':
                low_depth_dict[repid].append(sample_id)
    low_depth_json = json.dumps(low_depth_dict, indent=4)
    return low_depth_json


@click.command()
@click.option('--input-dir-eh', help='Input directory for ExpansionHunter VCFs')
@click.option(
    '--output-name-eh',
    help='Output file name for filter extractor output file eg eh.json',
)
def main(
    input_dir_eh,
    output_name_eh,
):
    """
    Hail batch job
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
    j.storage('20Gi')

    eh_json = j.call(eh_filter_extractor, input_dir_eh)

    b.write_output(eh_json.as_str(), output_path(output_name_eh, 'analysis'))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
