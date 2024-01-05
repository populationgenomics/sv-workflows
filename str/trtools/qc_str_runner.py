#!/usr/bin/env python3
# pylint: disable=duplicate-code
"""
This script runs qcSTR() from TRTools package on a single/merged STR vcf file and outputs various QC graphs

For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'tester' qcSTR_runner.py --caller=eh --file-path=gs://cpg-tob-wgs-test/hoptan-str/mergeSTR/mergeSTR_2_samples_gangstr.vcf
--refbias-binsize=5 --refbias-metric=mean --refbias-mingts=100 --refbias-xrange-min=0 --rebias-xrange-max=10000

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click
"""

import click

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path

config = get_config()

TRTOOLS_IMAGE = config['images']['trtools']


# inputs:
# file-path
@click.option('--file-path', help='gs://...')
# caller
@click.option(
    '--caller',
    help='gangstr or eh',
    type=click.Choice(['eh', 'gangstr'], case_sensitive=True),
)
# reference bias plot options
@click.option(
    '--refbias-binsize',
    help=' Sets the binsize (in bp) used to bin x-axis values, which give the reference TR length',
    default=5,
    type=int,
)
@click.option(
    '--refbias-metric',
    help='Determines which metric to use to summarize the reference bias in each bin.',
    type=click.Choice(['mean', 'median'], case_sensitive=True),
    default='mean',
)
@click.option(
    '--refbias-mingts',
    help='Exclude points computed using fewer than this many genotypes. This option is meant to avoid plotting outlier points driven by bins with small numbers of TRs with that reference length.',
    type=int,
    default=100,
)
@click.option(
    '--refbias-xrange-min',
    help='Exclude points corresponding to TRs with reference length less than this value.',
    type=int,
    default=0,
)
@click.option(
    '--refbias-xrange-max',
    help='Exclude points corresponding to TRs with reference length greater than this value.',
    type=int,
    default=10000,
)
@click.command()
def main(
    file_path,
    caller,
    refbias_binsize,
    refbias_metric,
    refbias_mingts,
    refbias_xrange_min,
    refbias_xrange_max,
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()
    vcf_input = b.read_input(file_path)
    trtools_job = b.new_job(name=f'qcSTR {caller}')

    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.storage('30G')
    trtools_job.cpu(8)

    if caller == 'eh':
        trtools_job.declare_resource_group(
            ofile={
                'sample-callnum.pdf': '{root}-sample-callnum.pdf',
                'chrom-callnum.pdf': '{root}-chrom-callnum.pdf',
                'diffref-histogram.pdf': '{root}-diffref-histogram.pdf',
                'diffref-bias.pdf': '{root}-diffref-bias.pdf'
                # EH does not have quality plots
            }
        )
        trtools_job.command(
            f"""
        set -ex;
        qcSTR --vcf {vcf_input} --vcftype {caller} --refbias-binsize {refbias_binsize} --refbias-metric {refbias_metric} --refbias-mingts {refbias_mingts} --refbias-xrange-min {refbias_xrange_min} --refbias-xrange-max {refbias_xrange_max} --out {trtools_job.ofile}

        """
        )

    elif caller == 'gangstr':
        trtools_job.declare_resource_group(
            ofile={
                'sample-callnum.pdf': '{root}-sample-callnum.pdf',
                'chrom-callnum.pdf': '{root}-chrom-callnum.pdf',
                'diffref-histogram.pdf': '{root}-diffref-histogram.pdf',
                'diffref-bias.pdf': '{root}-diffref-bias.pdf',
                'quality-per-locus.pdf': '{root}-quality-per-locus.pdf',
                'quality-sample-stratified.pdf': '{root}-quality-sample-stratified.pdf',
                'quality-per-sample.pdf': '{root}-quality-per-sample.pdf',
            }
        )
        trtools_job.command(
            f"""
        set -ex;
        qcSTR --vcf {vcf_input} --vcftype {caller} --quality per-locus --quality sample-stratified --quality per-sample --refbias-binsize {refbias_binsize} --refbias-metric {refbias_metric} --refbias-mingts {refbias_mingts} --refbias-xrange-min {refbias_xrange_min} --refbias-xrange-max {refbias_xrange_max} --out {trtools_job.ofile}

        """
        )
    else:
        raise ValueError('Invalid caller')

    output_path_vcf = output_path(f'qCSTR_samples_{caller}', 'analysis')
    b.write_output(trtools_job.ofile, output_path_vcf)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
