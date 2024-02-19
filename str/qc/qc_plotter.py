#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script makes QC plots

analysis-runner access-level "test" dataset "bioheart" description "QC plotter" output-dir "str/polymorphic_run/QC_plots" qc_plotter.py

"""
import hail as hl

from cpg_utils.config import get_config

from cpg_utils.hail_batch import init_batch, output_path
from bokeh.plotting import output_file, save

config = get_config()

def main():
    init_batch(worker_memory='highmem')
    mt = hl.read_matrix_table('gs://cpg-bioheart-test/str/polymorphic_run/mt_joiner/n_2045_annotated.mt')
    print(f'MT dimensions: {mt.count()}')

    non_mode_alleles_mt = mt.select_rows(
    allele_minus_mode = hl.agg.collect(mt.allele_1_minus_mode)
        .extend(hl.agg.collect(mt.allele_2_minus_mode))
    ).rows()
    non_mode_alleles_mt = non_mode_alleles_mt.explode('allele_minus_mode', name='alleles_minus_mode')
    p = hl.plot.histogram(non_mode_alleles_mt.alleles_minus_mode,legend= "Allele sizes minus mode", range=(-50, 100))

    output_file('local_plot.html')
    save(p)
    gcs_path = output_path('allele_size_range.html', 'analysis')
    hl.hadoop_copy('local_plot.html', gcs_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter