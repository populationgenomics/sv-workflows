#!/usr/bin/env python3
# pylint: disable=duplicate-code
"""
This script runs compareSTR() from TRTools package to compare calls between 2 STR VCFs.

For example: 
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' compareSTR_runner.py --file-path-1=gs://cpg-tob-wgs-test/hoptan-str/mergeSTR/mergeSTR_2_samples_gangstr.vcf
--file-path-2=gs://cpg-tob-wgs-test/hoptan-str/mergeSTR/mergeSTR_2_samples_eh.vcf --caller-1=gangstr --caller-2=eh

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""
import os
import click
import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path

config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']
BCFTOOLS_IMAGE = config['images']['bcftools']


# inputs:
# file-path-1
@click.option('--file-path-1', help='gs://...to VCF 1')
# file-path-2
@click.option('--file-path-2', help='gs://... to VCF 2')
# caller-1
@click.option('--caller-1', help='gangstr or eh')
# caller-2
@click.option('--caller-2', help='gangstr or eh')
@click.command()
def main(
    file_path_1, file_path_2, caller_1, caller_2
):  # pylint: disable=missing-function-docstring

    # Initializing Batch
    backend = hb.ServiceBackend(
        billing_project=get_config()['hail']['billing_project'],
        remote_tmpdir=remote_tmpdir(),
    )
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))
    vcf_input_1 = b.read_input(file_path_1)
    vcf_input_2 = b.read_input(file_path_2)

    # BCFTools is needed to sort, zip, and index files before using them as input in TRTools
    bcftools_job = b.new_job('Files prep')
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.storage('20G')
    bcftools_job.cpu(8)

    # declare resource groups, including extensions
    bcftools_job.declare_resource_group(
        vcf_1={'vcf.gz': '{root}.vcf_1.vcf.gz', 'vcf.gz.tbi': '{root}.vcf_1.vcf.gz.tbi'}
    )
    bcftools_job.declare_resource_group(
        vcf_2={'vcf.gz': '{root}.vcf_2.vcf.gz', 'vcf.gz.tbi': '{root}.vcf_2.vcf.gz.tbi'}
    )

    bcftools_job.command(
        f"""
    set -ex;
    
    echo "compressing {vcf_input_1}";
    bcftools sort {vcf_input_1} | bgzip -c >{bcftools_job.vcf_1['vcf.gz']};
    
    echo "indexing {bcftools_job.vcf_1['vcf.gz']}";
    tabix -p vcf {bcftools_job.vcf_1['vcf.gz']};
    
    echo "compressing {vcf_input_2}"; 
    bcftools sort {vcf_input_2} | bgzip -c >{bcftools_job.vcf_2['vcf.gz']};
    
    echo "indexing {bcftools_job.vcf_2['vcf.gz']}";
    tabix -p vcf {bcftools_job.vcf_2['vcf.gz']};
    """
    )

    trtools_job = b.new_job(name=f'compareSTR')
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('20G')
    trtools_job.cpu(8)

    trtools_job.declare_resource_group(
        ofile={
            'overall.tab': '{root}-overall.tab',
            'bubble-periodALL.pdf': '{root}-bubble-periodALL.pdf',
            'locuscompare.tab': '{root}-locuscompare.tab',
            'locuscompare.pdf': '{root}-locuscompare.pdf',
            'samplecompare.tab': '{root}-samplecompare.tab',
            'samplecompare.pdf': '{root}-samplecompare.pdf',
        }
    )

    trtools_job.command(
        f"""
    set -ex;
    compareSTR --vcf1 {bcftools_job.vcf_1['vcf.gz']} --vcf2 {bcftools_job.vcf_2['vcf.gz']} --vcftype1 {caller_1} --vcftype2 {caller_2} --out {trtools_job.ofile}
    
    """
    )
    output_path_vcf = output_path(f'compareSTR_samples_{caller_1}_{caller_2}')
    b.write_output(trtools_job.ofile, output_path_vcf)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
