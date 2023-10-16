#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script is step 3 of 4 for running associaTR.
It aims to:
- apply locus filters to mergedSTR VCF
- bgzip and tabix the filtered VCF for input into associatr
- remove "CPG" prefix from CPG IDs in the VCF
*only has to be run once per VCF file

 analysis-runner --dataset "tob-wgs" \
    --description "Run dumpSTR" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
    associatr_runner_part3_2.py --file-path=gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/dumpSTR/filtered_mergeSTR_results.filtered_vcf


"""
import click

from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch
from cpg_utils import to_path

from cpg_utils.hail_batch import output_path

config = get_config()
TRTOOLS_IMAGE = config['images']['trtools']
BCFTOOLS_IMAGE = config['images']['bcftools']

def file_parser(file_path, ofile_path):
    output_file = "relabelled.vcf"
    with to_path(file_path).open('r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            modified_line = line.replace("CPG", "")
            outfile.write(modified_line)


# inputs:
# file-path
@click.option('--file-path', help='gs://... to the output of mergedSTR')
@click.command()
def main(file_path):

    b = get_batch()
    file_parser_job = b.new_python_job(name = f'file_parser')
    file_parser_job.image(config['workflow']['driver_image'])
    file_parser_job.call(file_parser, file_path, file_parser_job.ofile)



    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter







