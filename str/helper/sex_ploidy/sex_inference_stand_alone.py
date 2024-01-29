#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script runs annotate_sex()

 analysis-runner --dataset "bioheart" \
    --description "standalone annotate-sex method" \
    --access-level "test" \
    --output-dir "qc-stand-alone/annotate-sex" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
    sex_inference_stand_alone.py --vds-path=gs://cpg-bioheart-test/vds/5-0.vds

"""

import logging

import hail as hl
from cpg_utils.hail_batch import reference_path, genome_build, output_path, init_batch
from gnomad.sample_qc.pipeline import annotate_sex
from hail.vds.variant_dataset import VariantDataset
from cpg_workflows.utils import can_reuse


import click


def generate_sex_coverage_mt(
    vds: hl.vds.VariantDataset,
    calling_intervals: hl.Table,
):
    """
    generate the pre-computed MT of coverage intervals, to supply to annotate_sex
    Args:
        vds ():
        calling_intervals ():

    Returns:

    """
    key_dtype = calling_intervals.key.dtype
    if (
        len(key_dtype) != 1
        or not isinstance(calling_intervals.key[0].dtype, hl.tinterval)
        or calling_intervals.key[0].dtype.point_type != vds.reference_data.locus.dtype
    ):
        raise ValueError(
            f'"impute_sex_chromosome_ploidy": expect calling_intervals to be list of intervals or '
            f'table with single key of type interval<locus>, found table with key: {key_dtype}'
        )

    rg = vds.reference_data.locus.dtype.reference_genome
    par_boundaries = []
    for par_interval in rg.par:
        par_boundaries.append(par_interval.start)
        par_boundaries.append(par_interval.end)

    # segment on PAR interval boundaries
    calling_intervals = hl.segment_intervals(calling_intervals, par_boundaries)

    # remove intervals overlapping PAR
    calling_intervals = calling_intervals.filter(
        hl.all(lambda x: ~x.overlaps(calling_intervals.interval), hl.literal(rg.par))
    )

    # checkpoint for efficient multiple downstream usages
    checkpoint_path = output_path(f'sample_qc/calling_intervals_partial.ht', 'tmp')
    if can_reuse(checkpoint_path, overwrite=True):
        calling_intervals = hl.read_matrix_table(str(checkpoint_path))
    else:
        calling_intervals = calling_intervals.checkpoint(str(checkpoint_path))

    interval = calling_intervals.key[0]
    (any_bad_intervals, chrs_represented) = calling_intervals.aggregate(
        (
            hl.agg.any(interval.start.contig != interval.end.contig),
            hl.agg.collect_as_set(interval.start.contig),
        )
    )
    if any_bad_intervals:
        raise ValueError(
            "'impute_sex_chromosome_ploidy' does not support calling intervals that span chromosome boundaries"
        )

    if len(rg.x_contigs) != 1:
        raise NotImplementedError(
            f"reference genome {rg.name!r} has multiple X contigs, this is not supported in 'impute_sex_chromosome_ploidy'"
        )
    if len(rg.y_contigs) != 1:
        raise NotImplementedError(
            f"reference genome {rg.name!r} has multiple Y contigs, this is not supported in 'impute_sex_chromosome_ploidy'"
        )

    kept_contig_filter = hl.array(chrs_represented).map(
        lambda x: hl.parse_locus_interval(x, reference_genome=rg)
    )
    vds = VariantDataset(
        hl.filter_intervals(vds.reference_data, kept_contig_filter),
        hl.filter_intervals(vds.variant_data, kept_contig_filter),
    )

    mt = vds.variant_data
    calling_intervals = calling_intervals.annotate(interval_dup=interval)
    mt = mt.annotate_rows(interval=calling_intervals[mt.locus].interval_dup)
    mt = mt.filter_rows(hl.is_defined(mt.interval))
    coverage = mt.select_entries(sum_dp=mt.DP, interval_size=hl.is_defined(mt.DP))
    coverage = coverage.key_rows_by('locus')

    # checkpoint the coverage mt prior to returning
    checkpoint_path = output_path('sample_qc/sex_coverage_precomputed.ht', 'tmp')

    if can_reuse(checkpoint_path, overwrite=True):
        coverage = hl.read_matrix_table(str(checkpoint_path))
    else:
        coverage = coverage.checkpoint(str(checkpoint_path))

    return coverage


def impute_sex(
    vds_path: str,
) -> hl.Table:
    """
    Impute sex based on coverage.
    """

    init_batch()

    # run()
    vds = hl.vds.read_vds(str(vds_path))

    # Remove centromeres and telomeres:
    tel_cent_ht = hl.read_table(str(reference_path('gnomad/tel_and_cent_ht')))
    if tel_cent_ht.count() > 0:
        vds = hl.vds.filter_intervals(vds, tel_cent_ht, keep=False)

    # Load calling intervals
    calling_intervals_path = reference_path(f'broad/genome_calling_interval_lists')
    calling_intervals_ht = hl.import_locus_intervals(
        str(calling_intervals_path), reference_genome=genome_build()
    )


    logging.info('Calling intervals table:')
    calling_intervals_ht.describe()

    # Pre-filter here and setting `variants_filter_lcr` and `variants_filter_segdup`
    # below to `False` to avoid the function calling gnomAD's `resources` module:
    for name in ['lcr_intervals_ht', 'seg_dup_intervals_ht']:
        interval_table = hl.read_table(str(reference_path(f'gnomad/{name}')))
        if interval_table.count() > 0:
            # remove all rows where the locus falls within a defined interval
            tmp_variant_data = vds.variant_data.filter_rows(
                hl.is_defined(interval_table[vds.variant_data.locus]), keep=False
            )
            vds = VariantDataset(
                reference_data=vds.reference_data, variant_data=tmp_variant_data
            ).checkpoint(output_path(f'{name}.vds', 'tmp'))
            logging.info(f'count post {name} filter:{vds.variant_data.count()}')

    # generate sex coverage mt, checkpointed
    coverage_mt = generate_sex_coverage_mt(vds, calling_intervals_ht)

    # Infer sex (adds row fields: is_female, var_data_chr20_mean_dp, sex_karyotype)
    sex_ht = annotate_sex(
        vds,
        included_intervals=calling_intervals_ht,
        gt_expr='LGT',
        variants_only_x_ploidy=True,
        variants_only_y_ploidy=False,
        coverage_mt=coverage_mt,
        variants_filter_lcr=False,  # already filtered above
        variants_filter_segdup=False,  # already filtered above
        variants_filter_decoy=False,
    )
    logging.info('Sex table:')
    sex_ht.describe()
    sex_ht = sex_ht.transmute(
        impute_sex_stats=hl.struct(
            f_stat=sex_ht.f_stat,
            n_called=sex_ht.n_called,
            expected_homs=sex_ht.expected_homs,
            observed_homs=sex_ht.observed_homs,
        )
    )

    # output writing
    sex_ht.write(output_path(f'sex.ht', 'analysis'))


@click.option(
    '--vds-path',
    help='GCS file path to VDS',
    type=str,
)
@click.command()
def main(vds_path):
    impute_sex(vds_path)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
