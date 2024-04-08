#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This Hail Query script outputs a histogram plotting mean chromosome (default chrY) coverage per individual using DP stats for each locus on
the specified chromosome (VDS file). Plots are stratified by inferred karyotype (X, XX, XY), provided with a supplementary sex-sample-mapping-file.

Input: VDS file path, sex-sample-mapping file path, and chromosome (default chrY)

 analysis-runner --dataset "bioheart" \
    --description "mean chromsome coverage plotter" \
    --access-level "test" \
    --output-dir "hoptan-str/mean_coverage_plot" \
    mean_coverage_plotter.py --vds-file-path=gs://cpg-bioheart-test/vds/5-0.vds --sex-sample-mapping-path=gs://cpg-bioheart-test/str/tester_sex_mapper.csv \
    --chromosome=chrY
"""
import csv

import click
from bokeh.plotting import output_file, save

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def coverage_plotter(vds_path, sex_sample_mapping, chromosome, xy_ylim):
    init_batch()
    vds = hl.vds.read_vds(vds_path)

    # extract variant matrix table from VDS
    mt_variant = vds.variant_data

    # print number of variants in VDS
    print(f'Dimensions of entire VDS: {mt_variant.count()}')

    # filter for chromosome:
    filtered_mt = mt_variant.filter_rows(mt_variant.locus.contig == chromosome)

    # separate samples by inferred karyotypic sex
    xx_samples = []
    xy_samples = []
    x_samples = []

    with to_path(sex_sample_mapping).open() as csv_file:
        csv_reader = csv.reader(csv_file)
        next(csv_reader)  # Skip the header row if it exists

        for row in csv_reader:
            sample_id = row[0]
            sex = row[2]

            if sex == 'XX':
                xx_samples.append(sample_id)
            elif sex == 'XY':
                xy_samples.append(sample_id)
            elif sex == 'X':
                x_samples.append(sample_id)

    # classify as 'set' type for Hail
    xx_samples = hl.literal(set(xx_samples))
    xy_samples = hl.literal(set(xy_samples))
    x_samples = hl.literal(set(x_samples))

    # filter for samples of interest
    filtered_mt_xx = filtered_mt.filter_cols(xx_samples.contains(filtered_mt.s))

    # calculate mean coverage per sample
    filtered_mt_xx = filtered_mt_xx.annotate_cols(
        dp_mean_cols=hl.agg.sum(filtered_mt_xx.DP) / filtered_mt_xx.count_rows(),
    )

    # Create a histogram using Hail's plotting functions
    p_xx = hl.plot.histogram(
        filtered_mt_xx.dp_mean_cols,
        legend=f'Mean DP per individual (XX) for {chromosome} loci',
    )
    output_file('local_plot_xx.html')
    save(p_xx)
    hl.hadoop_copy(
        'local_plot_xx.html',
        output_path(f'mean_coverage_xx_{chromosome}.html', 'analysis'),
    )

    # repeat for XY:
    filtered_mt_xy = filtered_mt.filter_cols(xy_samples.contains(filtered_mt.s))
    filtered_mt_xy = filtered_mt_xy.annotate_cols(
        dp_mean_cols=hl.agg.sum(filtered_mt_xy.DP) / filtered_mt_xy.count_rows(),
    )
    p_xy = hl.plot.histogram(
        filtered_mt_xy.dp_mean_cols,
        range=(0, xy_ylim),
        legend=f'Mean DP per individual (XY) for {chromosome} loci',
    )
    output_file('local_plot_xy.html')
    save(p_xy)
    hl.hadoop_copy(
        'local_plot_xy.html',
        output_path(f'mean_coverage_xy_{chromosome}.html', 'analysis'),
    )

    # repeat for X:
    filtered_mt_x = filtered_mt.filter_cols(x_samples.contains(filtered_mt.s))
    filtered_mt_x = filtered_mt_x.annotate_cols(dp_mean_cols=hl.agg.sum(filtered_mt_x.DP) / filtered_mt_x.count_rows())
    p_x = hl.plot.histogram(
        filtered_mt_x.dp_mean_cols,
        legend=f'Mean DP per individual (X) for {chromosome} loci',
    )
    output_file('local_plot_x.html')
    save(p_x)
    hl.hadoop_copy(
        'local_plot_x.html',
        output_path(f'mean_coverage_x_{chromosome}.html', 'analysis'),
    )


@click.option(
    '--vds-file-path',
    help='GCS file path to VDS file',
    type=str,
)
@click.option('--job-storage', help='Storage of the Hail batch job eg 30G', default='20G')
@click.option('--job-memory', help='Memory of the Hail batch job eg 64G', default='8G')
@click.option('--sex-sample-mapping-path', help='GCS file path to sex sample mapping file')
@click.option('--chromosome', help='Chromosome to plot mean coverage for')
@click.option('--xy-ylim', help='Y-axis limit for XY plot', default=100)
@click.command()
def main(vds_file_path, job_storage, job_memory, sex_sample_mapping_path, chromosome, xy_ylim):
    # Initialise batch
    b = get_batch()

    # Initialise job
    j = b.new_python_job(name=f'Mean coverage plotter job for {chromosome}')
    j.memory(job_memory)
    j.storage(job_storage)

    j.call(coverage_plotter, vds_file_path, sex_sample_mapping_path, chromosome, xy_ylim)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
