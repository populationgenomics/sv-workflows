#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
 analysis-runner --dataset "tob-wgs" \
    --description "Run dumpSTR" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
    associatr_runner_part3_2.py --file-path=gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/dumpSTR/filtered_mergeSTR_results.vcf
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
    with to_path(file_path).open('r') as infile, open(ofile_path, 'w') as outfile:
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

    bcftools_job = b.new_job(name = f'bgzip and tabix the dumpSTR output VCF')
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.storage('20G')
    bcftools_job.cpu(4)

    bcftools_job.declare_resource_group(
        vcf_output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    bcftools_job.command(
        f"""
    set -ex;
    echo "Compressing";
    bcftools sort {file_parser_job.ofile} | bgzip -c > {bcftools_job.vcf_output['vcf.gz']};

    echo "indexing {bcftools_job.vcf_output['vcf.gz']}";
    tabix -p vcf {bcftools_job.vcf_output['vcf.gz']};
"""
    )
    b.write_output(bcftools_job.vcf_output, output_path(f'input_files/dumpSTR/dumpSTR_filtered'))


    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter







