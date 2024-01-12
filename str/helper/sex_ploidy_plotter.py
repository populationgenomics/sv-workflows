#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This Hail Query script outputs a scatterplot plotting chrX and chrY ploidy, and labelling the points by inferred karyotypic sex from the SampleQC Hail Table.

 analysis-runner --dataset "bioheart" \
    --description "karyotypic sex extractor" \
    --access-level "test" \
    --output-dir "hoptan-str/sex_ploidy_plot" \
    sex_ploidy_plotter.py --file-path=gs://cpg-bioheart-test/large_cohort/1-0/sample_qc.ht/

"""

import hail as hl
import click

from cpg_utils.hail_batch import output_path, init_batch, get_batch
from bokeh.plotting import figure, output_file, save


def sex_ploidy_plotter(file_path):
    init_batch()
    sample_qc_table = hl.read_table(file_path)

    # Create a scatterplot matrix using Hail's plotting functions
    p = hl.plot.scatter(
        x=sample_qc_table.chrY_ploidy,
        y=sample_qc_table.chrX_ploidy,
        label=sample_qc_table.sex_karyotype,
        title='Scatterplot of chrY_ploidy vs chrX_ploidy',
        xlabel='chrY_ploidy',
        ylabel='chrX_ploidy'
    )

    # Save the plot to a local file, then hadoop_copy to copy to GCS bucket
    return save(p)
    #hl.hadoop_copy('local_plot.html', "gs://cpg-tob-wgs-test/hoptan-str/sex_ploidy_plot/sex_ploidy_plot.html")


@click.option(
    '--file-path',
    help='GCS file path to SampleQC Hail Table',
    type=str,
)
@click.command()
def main(file_path):

    b = get_batch()
    j = b.new_python_job()
    gcs_output_path = output_path(f'sex_ploidy_plot.html', 'analysis')
    result = j.call(sex_ploidy_plotter, file_path)
    b.write_output(result.as_str(), gcs_output_path)

    b.run(wait=False)
if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
