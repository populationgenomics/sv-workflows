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
    #print(f'MT dimensions: {mt.count()}')


    # Filter out monomorphic loci
    mt = mt.filter_rows(mt.num_alleles>1)
    #print(f'MT dimensions after subsetting to loci with more than 1 allele: {mt.count()}')


    # Alleles minus mode histogram
    alleles_minus_mode_ht = mt.select_rows(
    allele_minus_mode = hl.agg.collect(mt.allele_1_minus_mode)
        .extend(hl.agg.collect(mt.allele_2_minus_mode))
    ).rows()
    alleles_minus_mode_ht = alleles_minus_mode_ht.explode('allele_minus_mode', name='alleles_minus_mode')
    p = hl.plot.histogram(alleles_minus_mode_ht.alleles_minus_mode,legend= "Allele sizes minus mode", range=(-50, 100))

    output_file('local_plot.html')
    save(p)
    gcs_path = output_path('alleles_minus_mode/allele_size_range_50_to_100.html', 'analysis')
    hl.hadoop_copy('local_plot.html', gcs_path)

    q = hl.plot.histogram(alleles_minus_mode_ht.alleles_minus_mode,legend= "Allele sizes minus mode")
    output_file('local_plot_1.html')
    save(q)
    gcs_path_1 = output_path('alleles_minus_mode/allele_size_range.html', 'analysis')
    hl.hadoop_copy('local_plot_1.html', gcs_path_1)

    r = hl.plot.histogram(alleles_minus_mode_ht.alleles_minus_mode,legend= "Allele sizes minus mode", range=(-20, 20))
    output_file('local_plot_2.html')
    save(r)
    gcs_path_2 = output_path('alleles_minus_mode/allele_size_range_20_20.html', 'analysis')
    hl.hadoop_copy('local_plot_2.html', gcs_path_2)

    # Alleles histogram (num. repeats)
    alleles_ht = mt.select_rows(
    alleles = hl.agg.collect(mt.allele_1_rep_length)
        .extend(hl.agg.collect(mt.allele_2_rep_length))
    ).rows()
    alleles_ht = alleles_ht.explode('alleles', name='alleles')
    s = hl.plot.histogram(alleles_ht.alleles,legend= "Allele sizes (num. repeats)")
    output_file('local_plot_3.html')
    save(s)
    gcs_path_3 = output_path('alleles_num_repeats/alleles_num_repeats.html', 'analysis')
    hl.hadoop_copy('local_plot_3.html', gcs_path_3)

    t = hl.plot.histogram(alleles_ht.alleles,legend= "Allele sizes (num. repeats)", range=(0, 50))
    output_file('local_plot_4.html')
    save(t)
    gcs_path_4 = output_path('alleles_num_repeats/alleles_num_repeats_0_50.html', 'analysis')
    hl.hadoop_copy('local_plot_4.html', gcs_path_4)

    # Alleles histogram (bp length)
    alleles_bp_ht = mt.select_rows(
    alleles_bp = hl.agg.collect(mt.allele_1_bp_length)
        .extend(hl.agg.collect(mt.allele_2_bp_length)
    )
    ).rows()
    alleles_bp_ht = alleles_bp_ht.explode('alleles_bp', name='alleles_bp')
    u = hl.plot.histogram(alleles_bp_ht.alleles_bp,legend= "Allele sizes (bp length)")
    output_file('local_plot_5.html')
    save(u)
    gcs_path_5 = output_path('alleles_bp_length/alleles_bp_length.html', 'analysis')
    hl.hadoop_copy('local_plot_5.html', gcs_path_5)

    v = hl.plot.histogram(alleles_bp_ht.alleles_bp,legend= "Allele sizes (bp length)", range=(0, 500))
    output_file('local_plot_6.html')
    save(v)
    gcs_path_6 = output_path('alleles_bp_length/alleles_bp_length_0_500.html', 'analysis')
    hl.hadoop_copy('local_plot_6.html', gcs_path_6)




    """
    #MT dimensions
    print(f'Original MT dimensions: {mt.count()}')
    mt = mt.filter_entries(mt.num_alleles > 1)
    print(f'MT dimensions after subsetting to loci with more than 1 allele: {mt.count()}')
    mt = mt.filter_entries((mt.allele_1_minus_mode> -21) & (mt.allele_1_minus_mode<21) & (mt.allele_2_minus_mode>-21) & (mt.allele_2_minus_mode<21))
    print(f'MT dimensions after subsetting to [-20,20] alleles rel. to mode: {mt.count()}')
    """



if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter