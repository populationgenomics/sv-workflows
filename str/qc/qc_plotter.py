#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script makes QC plots

analysis-runner --access-level "test" --dataset "bioheart" --description "QC plotter" --output-dir "str/polymorphic_run_n2045/QC" qc_plotter.py \
--mt-path=gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt

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


    # Filter out monomorphic loci,locus-level call rate 0.9 threshold, obs_het >= 0.00995
    #mt = mt.filter_rows(
    #    (mt.num_alleles > 1) & (mt.variant_qc.call_rate >= 0.9) & (mt.obs_het >= 0.00995) & (mt.binom_hwep >= 0.000001),
    #)
    #print(f'MT dimensions after locus-level filters: {mt.count()}')
    # Filter out chrX
    #mt = mt.filter_rows(mt.locus.contig != 'chrX')
    #print(f'MT dimensions after filtering out chrX: {mt.count()}')
    #potato = mt.filter_entries((mt.allele_1_minus_mode> -21) & (mt.allele_1_minus_mode<21) & (mt.allele_2_minus_mode>-21) & (mt.allele_2_minus_mode<21))
    #print(f' MT cap [-20,20] rel. to mode: {potato.entries().count()}')
    #mt.rows().export('gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/annotated_mt/v2/str_annotated_rows.tsv.bgz')

    alleles_minus_mode_ht = mt.select_rows(
    allele_minus_mode = hl.agg.collect(mt.allele_1_minus_mode)
        .extend(hl.agg.collect(mt.allele_2_minus_mode))
    ).rows()
    alleles_minus_mode_ht = alleles_minus_mode_ht.explode('allele_minus_mode', name='alleles_minus_mode')
    #alleles_minus_mode_ht.export('gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/annotated_mt/v2/alleles_minus_mode_ht.tsv.bgz')
    #pq = hl.plot.histogram(alleles_minus_mode_ht.alleles_minus_mode,legend= "Allele relative to mode allele", range = (-20, 20))
    #output_file('local_plot_pq.html')
    #save(pq)
    #gcs_path_pq = output_path('alleles_minus_mode/range_20_20', 'analysis')
    #hl.hadoop_copy('local_plot_pq.html', gcs_path_pq)

    # calculate proportion of alleles that are within 20bp of the mode allele
    alleles_minus_mode_ht = alleles_minus_mode_ht.annotate(
    within_20bp = (alleles_minus_mode_ht.alleles_minus_mode >= -20) & (alleles_minus_mode_ht.alleles_minus_mode <= 20),
    within_10bp = (alleles_minus_mode_ht.alleles_minus_mode >= -10) & (alleles_minus_mode_ht.alleles_minus_mode <= 10),
)

    within_20bp_proportion = alleles_minus_mode_ht.aggregate(hl.agg.fraction(alleles_minus_mode_ht.within_20bp))
    print(f'Proportion of alleles within 20bp of the mode allele: {within_20bp_proportion}')

    within_10bp_proportion = alleles_minus_mode_ht.aggregate(hl.agg.fraction(alleles_minus_mode_ht.within_10bp))
    print(f'Proportion of alleles within 10bp of the mode allele: {within_10bp_proportion}')


    #chr = 'chr14'
    #position = 42532368

    #mt = mt.filter_rows((mt.locus.contig == chr) & (mt.locus.position == position))

    #table = hl.import_table('gs://cpg-bioheart-test/str/sample-sex-mapping/sample_karyotype_sex_mapping.csv', delimiter=',', impute=True)
    # Define a function to determine the cohort based on the 'id' column
    #def get_cohort(id):
    #    return hl.if_else(hl.str(id).startswith('CT'), 'bioheart', 'tob')

    # Add a new column 'cohort' based on the 'id' column using the defined function
    #table = table.annotate(cohort=get_cohort(table.external_id))
    #table = table.key_by('s')
    #mt = mt.annotate_cols(cohort = table[mt.s].cohort)


    #t = mt.annotate_entries(allele_1_minus_ref = mt.allele_1_rep_length- mt.info.REF)
    #t = mt.annotate_entries(allele_2_minus_ref = mt.allele_2_rep_length-mt.info.REF)

    #dotato = mt.filter_entries((mt.allele_1_minus_ref> -21) & (mt.allele_1_minus_ref<21) & (mt.allele_2_minus_ref>-21) & (mt.allele_2_minus_ref<21))
    #print(f' MT cap [-20,20] rel. to ref: {dotato.entries().count()}')

    # Alleles minus ref histogram
    ##alleles_minus_ref_ht = mt.select_rows(
    #allele_minus_ref = hl.agg.collect(mt.allele_1_minus_mode)
    #    .extend(hl.agg.collect(mt.allele_2_minus_mode))
    #).rows()
    #alleles_minus_ref_ht = alleles_minus_ref_ht.explode('allele_minus_ref', name='alleles_minus_ref')
    #p = hl.plot.histogram(alleles_minus_ref_ht.alleles_minus_ref,legend= "Allele sizes minus mode (10to200)", range = (10,200))
    #output_file('local_plot_1.html')
    #save(p)
    #gcs_path_1 = output_path(f'alleles_minus_mode/10to200.html', 'analysis')
    #hl.hadoop_copy('local_plot_1.html', gcs_path_1)

    #for cohort in ['tob','bioheart']:
        #mt_cohort = mt.filter_cols(mt.cohort == cohort)

        # Alleles minus ref histogram


    """
    # Alleles minus mode histogram
    mt = mt.filter_rows(mt.locus == hl.locus('chr1', 105963350))
    alleles_minus_mode_ht = mt.select_rows(
    allele_minus_mode = hl.agg.collect(mt.allele_1_minus_mode)
        .extend(hl.agg.collect(mt.allele_2_minus_mode))
    ).rows()
    alleles_minus_mode_ht = alleles_minus_mode_ht.explode('allele_minus_mode', name='alleles_minus_mode')

    pq = hl.plot.histogram(alleles_minus_mode_ht.alleles_minus_mode,legend= "Allele sizes minus mode, chr1:105963350")
    output_file('local_plot_pq.html')
    save(pq)
    gcs_path_pq = output_path('alleles_minus_mode/allele_size_chr1_105963350.html', 'analysis')
    hl.hadoop_copy('local_plot_pq.html', gcs_path_pq)


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
    """
    """
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

    t = hl.plot.histogram(alleles_ht.alleles,legend= "Allele sizes (num. repeats)", range=(0, 100))
    output_file('local_plot_4.html')
    save(t)
    gcs_path_4 = output_path('alleles_num_repeats/alleles_num_repeats_0_100.html', 'analysis')
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
    ## CI over CN plots

    CI_over_CN_ht = mt.select_rows(
    CI_over_CN = hl.agg.collect(mt.allele_1_REPCI_over_CN)
        .extend(hl.agg.collect(mt.allele_2_REPCI_over_CN)
    )
    ).rows()
    CI_over_CN_ht = CI_over_CN_ht.explode('CI_over_CN', name='CI_over_CN')
    w = hl.plot.histogram(CI_over_CN_ht.CI_over_CN,legend= "CI/CN of all alleles")
    output_file('local_plot_7.html')
    save(w)
    gcs_path_7 = output_path('CI_over_CN/CI_over_CN.html', 'analysis')
    hl.hadoop_copy('local_plot_7.html', gcs_path_7)

    x = hl.plot.histogram(CI_over_CN_ht.CI_over_CN,legend= "CI/CN of all alleles", range=(0, 5))
    output_file('local_plot_8.html')
    save(x)
    gcs_path_8 = output_path('CI_over_CN/CI_over_CN_0_5.html', 'analysis')
    hl.hadoop_copy('local_plot_8.html', gcs_path_8)

    y = hl.plot.histogram(CI_over_CN_ht.CI_over_CN,legend= "CI/CN of all alleles", range=(0, 2))
    output_file('local_plot_9.html')
    save(y)
    gcs_path_9 = output_path('CI_over_CN/CI_over_CN_0_2.html', 'analysis')
    hl.hadoop_copy('local_plot_9.html', gcs_path_9)

    z = hl.plot.histogram(CI_over_CN_ht.CI_over_CN,legend= "CI/CN of all alleles", range=(0, 1))
    output_file('local_plot_10.html')
    save(z)
    gcs_path_10 = output_path('CI_over_CN/CI_over_CN_0_1.html', 'analysis')
    hl.hadoop_copy('local_plot_10.html', gcs_path_10)


    #MT dimensions
    print(f'Original MT dimensions: {mt.count()}')
    mt = mt.filter_entries(mt.num_alleles > 1)
    print(f'MT dimensions after subsetting to loci with more than 1 allele: {mt.count()}')
    mt = mt.filter_entries((mt.allele_1_minus_mode> -21) & (mt.allele_1_minus_mode<21) & (mt.allele_2_minus_mode>-21) & (mt.allele_2_minus_mode<21))
    print(f'MT dimensions after subsetting to [-20,20] alleles rel. to mode: {mt.count()}')
    """



if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter