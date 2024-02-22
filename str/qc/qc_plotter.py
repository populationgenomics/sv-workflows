#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script makes QC plots

analysis-runner --access-level "test" --dataset "bioheart" --description "QC plotter" --output-dir "str/polymorphic_run_n2045/QC" qc_plotter.py \
--mt-path=gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v1/str_annotated.mt

"""
import hail as hl
import click

from cpg_utils.config import get_config

from cpg_utils.hail_batch import init_batch, output_path
from bokeh.plotting import output_file, save

config = get_config()

@click.option('--mt-path', help='GCS Path to the input MT')
@click.command()
def main(mt_path):
    init_batch()

    mt = hl.read_matrix_table(mt_path)
    print(f'MT dimensions: {mt.count()}')

    # Filter out monomorphic loci
    mt = mt.filter_rows(mt.num_alleles>1)
    print(f'MT dimensions after subsetting to loci with more than 1 allele: {mt.count()}')

    # Locus level call rate
    locus_call_rate = hl.plot.histogram(mt.variant_qc.call_rate,legend= "Locus level call rate")
    output_file('locus_call_rate.html')
    save(locus_call_rate)
    gcs_path = output_path('locus_call_rate/locus_call_rate.html', 'analysis')
    hl.hadoop_copy('locus_call_rate.html', gcs_path)




    """
    non_mode_alleles_mt = mt.select_rows(
    allele_minus_mode = hl.agg.collect(mt.allele_1_minus_mode)
        .extend(hl.agg.collect(mt.allele_2_minus_mode))
    ).rows()
    non_mode_alleles_mt = non_mode_alleles_mt.explode('allele_minus_mode', name='alleles_minus_mode')
    p = hl.plot.histogram(non_mode_alleles_mt.alleles_minus_mode,legend= "Allele sizes minus mode", range=(-50, 100))

    output_file('local_plot.html')
    save(p)
    gcs_path = output_path('allele_size_range/allele_size_range_50_to_100.html', 'analysis')
    hl.hadoop_copy('local_plot.html', gcs_path)

    q = hl.plot.histogram(non_mode_alleles_mt.alleles_minus_mode,legend= "Allele sizes minus mode")
    output_file('local_plot_1.html')
    save(q)
    gcs_path_1 = output_path('allele_size_range/allele_size_range.html', 'analysis')
    hl.hadoop_copy('local_plot_1.html', gcs_path_1)

    r = hl.plot.histogram(non_mode_alleles_mt.alleles_minus_mode,legend= "Allele sizes minus mode", range=(-20, 20))
    output_file('local_plot_2.html')
    save(r)
    gcs_path_2 = output_path('allele_size_range/allele_size_range_20_20.html', 'analysis')
    hl.hadoop_copy('local_plot_2.html', gcs_path_2)


    #MT dimensions
    print(f'Original MT dimensions: {mt.count()}')
    mt = mt.filter_entries(mt.num_alleles > 1)
    print(f'MT dimensions after subsetting to loci with more than 1 allele: {mt.count()}')
    mt = mt.filter_entries((mt.allele_1_minus_mode> -21) & (mt.allele_1_minus_mode<21) & (mt.allele_2_minus_mode>-21) & (mt.allele_2_minus_mode<21))
    print(f'MT dimensions after subsetting to [-20,20] alleles rel. to mode: {mt.count()}')
    """



if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter