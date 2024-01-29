#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script runs annotate_sex()

 analysis-runner --dataset "bioheart" \
    --description "standalone annotate-sex method" \
    --access-level "test" \
    --output-dir "qc-stand-alone/annotate-sex" \
    sex_inference_stand_alone.py --input-dir=gs://cpg-bioheart-test/vds/5-0.vds

"""

import logging

import hail as hl
from cpg_utils.config import get_config
from cpg_utils.hail_batch import reference_path, genome_build, output_path, init_batch
from hail.vds.variant_dataset import VariantDataset

import click

def annotate_sex(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    is_sparse: bool = True,
    excluded_intervals: Optional[hl.Table] = None,
    included_intervals: Optional[hl.Table] = None,
    normalization_contig: str = "chr20",
    sites_ht: Optional[hl.Table] = None,
    aaf_expr: Optional[str] = None,
    gt_expr: str = "GT",
    f_stat_cutoff: float = 0.5,
    aaf_threshold: float = 0.001,
    variants_only_x_ploidy: bool = False,
    variants_only_y_ploidy: bool = False,
    variants_filter_lcr: bool = True,
    variants_filter_segdup: bool = True,
    variants_filter_decoy: bool = False,
    variants_snv_only: bool = False,
    coverage_mt: Optional[hl.MatrixTable] = None,
    compute_x_frac_variants_hom_alt: bool = False,
    compute_fstat: bool = True,
    infer_karyotype: bool = True,
    use_gaussian_mixture_model: bool = False,
) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    Return Table with the following fields:
        - s (str): Sample
        - `normalization_contig`_mean_dp (float32): Sample's mean coverage over the specified `normalization_contig`.
        - chrX_mean_dp (float32): Sample's mean coverage over chromosome X.
        - chrY_mean_dp (float32): Sample's mean coverage over chromosome Y.
        - chrX_ploidy (float32): Sample's imputed ploidy over chromosome X.
        - chrY_ploidy (float32): Sample's imputed ploidy over chromosome Y.

        If `compute_fstat`:
            - f_stat (float64): Sample f-stat. Calculated using hl.impute_sex.
            - n_called (int64): Number of variants with a genotype call. Calculated using hl.impute_sex.
            - expected_homs (float64): Expected number of homozygotes. Calculated using hl.impute_sex.
            - observed_homs (int64): Observed number of homozygotes. Calculated using hl.impute_sex.

        If `infer_karyotype`:
            - X_karyotype (str): Sample's chromosome X karyotype.
            - Y_karyotype (str): Sample's chromosome Y karyotype.
            - sex_karyotype (str): Sample's sex karyotype.

    .. note::

            In order to infer sex karyotype (`infer_karyotype`=True), one of `compute_fstat` or
            `use_gaussian_mixture_model` must be set to True.

    :param mtds: Input MatrixTable or VariantDataset.
    :param is_sparse: Whether input MatrixTable is in sparse data format. Default is True.
    :param excluded_intervals: Optional table of intervals to exclude from the computation. This option is currently
        not implemented for imputing sex chromosome ploidy on a VDS.
    :param included_intervals: Optional table of intervals to use in the computation. REQUIRED for exomes.
    :param normalization_contig: Which chromosome to use to normalize sex chromosome coverage. Used in determining sex
        chromosome ploidies. Default is "chr20".
    :param sites_ht: Optional Table of sites and alternate allele frequencies for filtering the input MatrixTable prior to imputing sex.
    :param aaf_expr: Optional. Name of field in input MatrixTable with alternate allele frequency.
    :param gt_expr: Name of entry field storing the genotype. Default is 'GT'.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY
        samples are above cutoff. Default is 0.5.
    :param aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations. Default is 0.001.
    :param variants_only_x_ploidy: Whether to use depth of only variant data for the x ploidy estimation.
    :param variants_only_y_ploidy: Whether to use depth of only variant data for the y ploidy estimation.
    :param variants_filter_lcr: Whether to filter out variants in LCR regions for variants only ploidy estimation and
        fraction of homozygous alternate variants on chromosome X. Default is True.
    :param variants_filter_segdup: Whether to filter out variants in segdup regions for variants only ploidy estimation
        and fraction of homozygous alternate variants on chromosome X. Default is True.
    :param variants_filter_decoy: Whether to filter out variants in decoy regions for variants only ploidy estimation
        and fraction of homozygous alternate variants on chromosome X. Default is False. Note: this option doesn't
        exist for GRCh38.
    :param variants_snv_only: Whether to filter to only single nucleotide variants for variants only ploidy estimation
        and fraction of homozygous alternate variants on chromosome X. Default is False.
    :param coverage_mt: Optional precomputed coverage MatrixTable to use in reference based VDS ploidy estimation.
    :param compute_x_frac_variants_hom_alt: Whether to return an annotation for the fraction of homozygous alternate
        variants on chromosome X. Default is False.
    :param compute_fstat: Whether to compute f-stat. Default is True.
    :param infer_karyotype: Whether to infer sex karyotypes. Default is True.
    :param use_gaussian_mixture_model: Whether to use gaussian mixture model to split samples into 'XX' and 'XY'
        instead of f-stat. Default is False.
    :return: Table of samples and their imputed sex karyotypes.
    """
    logger.info("Imputing sex chromosome ploidies...")

    if infer_karyotype and not (compute_fstat or use_gaussian_mixture_model):
        raise ValueError(
            "In order to infer sex karyotype (infer_karyotype=True), one of"
            " 'compute_fstat' or 'use_gaussian_mixture_model' must be set to True!"
        )

    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        if excluded_intervals is not None:
            raise NotImplementedError(
                "The use of the parameter 'excluded_intervals' is currently not"
                " implemented for imputing sex chromosome ploidy on a VDS!"
            )
        if included_intervals is None:
            raise NotImplementedError(
                "The current implementation for imputing sex chromosome ploidy on a VDS"
                " requires a list of 'included_intervals'!"
            )
        mt = mtds.variant_data
    else:
        if not is_sparse:
            raise NotImplementedError(
                "Imputing sex ploidy does not exist yet for dense data."
            )
        mt = mtds

    # Determine the contigs that are needed for variant only and reference
    # block only sex ploidy imputation
    rg = get_reference_genome(mt.locus)
    if normalization_contig not in rg.contigs:
        raise ValueError(
            f"Normalization contig {normalization_contig} is not found in reference"
            f" genome {rg.name}!"
        )

    x_contigs = set(rg.x_contigs)
    y_contigs = set(rg.y_contigs)
    if variants_only_x_ploidy:
        var_keep_contigs = x_contigs | {normalization_contig}
        ref_keep_contigs = set()
    else:
        ref_keep_contigs = x_contigs | {normalization_contig}
        var_keep_contigs = set()
    if variants_only_y_ploidy:
        var_keep_contigs = {normalization_contig} | var_keep_contigs | y_contigs
    else:
        ref_keep_contigs = {normalization_contig} | ref_keep_contigs | y_contigs

    ref_keep_locus_intervals = [
        hl.parse_locus_interval(contig, reference_genome=rg.name)
        for contig in ref_keep_contigs
    ]
    var_keep_locus_intervals = [
        hl.parse_locus_interval(contig, reference_genome=rg.name)
        for contig in var_keep_contigs
    ]
    x_locus_intervals = [
        hl.parse_locus_interval(contig, reference_genome=rg.name)
        for contig in x_contigs
    ]

    if ref_keep_contigs:
        logger.info(
            "Imputing sex chromosome ploidy using only reference block depth"
            " information on the following contigs: %s",
            ref_keep_contigs,
        )
        if is_vds:
            if coverage_mt is not None:
                ploidy_ht = hl.vds.impute_sex_chr_ploidy_from_interval_coverage(
                    coverage_mt.filter_rows(
                        hl.is_defined(included_intervals[coverage_mt.row_key])
                        & hl.literal(ref_keep_contigs).contains(
                            coverage_mt.interval.start.contig
                        )
                    ),
                    normalization_contig=normalization_contig,
                )
            else:
                ploidy_ht = hl.vds.impute_sex_chromosome_ploidy(
                    hl.vds.filter_intervals(mtds, ref_keep_locus_intervals),
                    calling_intervals=included_intervals,
                    normalization_contig=normalization_contig,
                    use_variant_dataset=False,
                )
            ploidy_ht = ploidy_ht.rename(
                {
                    "x_ploidy": "chrX_ploidy",
                    "y_ploidy": "chrY_ploidy",
                    "x_mean_dp": "chrX_mean_dp",
                    "y_mean_dp": "chrY_mean_dp",
                }
            )
        else:
            ploidy_ht = impute_sex_ploidy(
                hl.filter_intervals(mt, ref_keep_locus_intervals),
                excluded_intervals,
                included_intervals,
                normalization_contig,
                use_only_variants=False,
            )
        if variants_only_x_ploidy:
            ploidy_ht = ploidy_ht.drop("chrX_ploidy", "chrX_mean_dp")
        if variants_only_y_ploidy:
            ploidy_ht = ploidy_ht.drop("chrY_ploidy", "chrY_mean_dp")

    add_globals = hl.struct()
    if compute_x_frac_variants_hom_alt or var_keep_contigs:
        logger.info(
            "Filtering variants for variant only sex chromosome ploidy imputation"
            " and/or computation of the fraction of homozygous alternate variants on"
            " chromosome X",
        )
        filtered_mt = hl.filter_intervals(
            mt, var_keep_locus_intervals + x_locus_intervals
        )
        if variants_filter_lcr or variants_filter_segdup or variants_filter_decoy:
            logger.info(
                "Filtering out variants in: %s",
                ("segmental duplications, " if variants_filter_segdup else "")
                + ("low confidence regions, " if variants_filter_lcr else "")
                + (" decoy regions" if variants_filter_decoy else ""),
            )
            filtered_mt = filter_low_conf_regions(
                filtered_mt,
                filter_lcr=variants_filter_lcr,
                filter_decoy=variants_filter_decoy,
                filter_segdup=variants_filter_segdup,
            )
        if variants_snv_only:
            logger.info("Filtering to SNVs")
            filtered_mt = filtered_mt.filter_rows(
                hl.is_snp(filtered_mt.alleles[0], filtered_mt.alleles[1])
            )

        add_globals = add_globals.annotate(
            variants_filter_lcr=variants_filter_lcr,
            variants_segdup=variants_filter_segdup,
            variants_filter_decoy=variants_filter_decoy,
            variants_snv_only=variants_snv_only,
        )

    if var_keep_contigs:
        logger.info(
            "Imputing sex chromosome ploidy using only variant depth information on the"
            " following contigs: %s",
            var_keep_contigs,
        )
        var_filtered_mt = hl.filter_intervals(filtered_mt, var_keep_locus_intervals)
        if is_vds:
            var_ploidy_ht = hl.vds.impute_sex_chromosome_ploidy(
                hl.vds.VariantDataset(mtds.reference_data, var_filtered_mt),
                calling_intervals=included_intervals,
                normalization_contig=normalization_contig,
                use_variant_dataset=True,
            )
            var_ploidy_ht = var_ploidy_ht.rename(
                {
                    "autosomal_mean_dp": f"var_data_{normalization_contig}_mean_dp",
                    "x_ploidy": "chrX_ploidy",
                    "y_ploidy": "chrY_ploidy",
                    "x_mean_dp": "chrX_mean_dp",
                    "y_mean_dp": "chrY_mean_dp",
                }
            )
        else:
            var_ploidy_ht = impute_sex_ploidy(
                var_filtered_mt,
                excluded_intervals,
                included_intervals,
                normalization_contig,
                use_only_variants=True,
            )
            var_ploidy_ht = var_ploidy_ht.rename(
                {
                    f"{normalization_contig}_mean_dp": (
                        f"var_data_{normalization_contig}_mean_dp"
                    )
                }
            )

        if ref_keep_contigs:
            ploidy_ht = var_ploidy_ht.annotate(**ploidy_ht[var_ploidy_ht.key])
        else:
            ploidy_ht = var_ploidy_ht

    ploidy_ht = ploidy_ht.annotate_globals(
        normalization_contig=normalization_contig,
        variants_only_x_ploidy=variants_only_x_ploidy,
        variants_only_y_ploidy=variants_only_y_ploidy,
        **add_globals,
    )

    if compute_x_frac_variants_hom_alt:
        logger.info(
            "Computing fraction of variants that are homozygous alternate on"
            " chromosome X"
        )
        filtered_mt = hl.filter_intervals(filtered_mt, x_locus_intervals)
        filtered_mt = filtered_mt.filter_rows(
            hl.is_defined(included_intervals[filtered_mt.locus])
        )
        filtered_mt = filtered_mt.annotate_entries(
            adj=get_adj_expr(
                filtered_mt.LGT, filtered_mt.GQ, filtered_mt.DP, filtered_mt.LAD
            )
        )
        frac_hom_alt_ht = filtered_mt.select_cols(
            chrx_frac_hom_alt=hl.agg.count_where(filtered_mt.LGT.is_hom_var())
            / hl.agg.count_where(hl.is_defined(filtered_mt.LGT)),
            chrx_frac_hom_alt_adj=hl.agg.filter(
                filtered_mt.adj,
                hl.agg.count_where(filtered_mt.LGT.is_hom_var())
                / hl.agg.count_where(hl.is_defined(filtered_mt.LGT)),
            ),
        ).cols()
        ploidy_ht = ploidy_ht.annotate(**frac_hom_alt_ht[ploidy_ht.key])

    if compute_fstat:
        logger.info("Filtering mt to biallelic SNPs in X contigs: %s", x_contigs)
        if "was_split" in list(mt.row):
            mt = mt.filter_rows(
                (~mt.was_split) & hl.is_snp(mt.alleles[0], mt.alleles[1])
            )
        else:
            mt = mt.filter_rows(
                (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
            )

        mt = hl.filter_intervals(mt, x_locus_intervals)
        if sites_ht is not None:
            if aaf_expr is None:
                logger.warning(
                    "sites_ht was provided, but aaf_expr is missing. Assuming name of"
                    " field with alternate allele frequency is 'AF'."
                )
                aaf_expr = "AF"
            logger.info("Filtering to provided sites")
            mt = mt.annotate_rows(**sites_ht[mt.row_key])
            mt = mt.filter_rows(hl.is_defined(mt[aaf_expr]))

        logger.info("Calculating inbreeding coefficient on chrX")
        sex_ht = hl.impute_sex(
            mt[gt_expr],
            aaf_threshold=aaf_threshold,
            male_threshold=f_stat_cutoff,
            female_threshold=f_stat_cutoff,
            aaf=aaf_expr,
        )

        logger.info("Annotating sex chromosome ploidy HT with impute_sex HT")
        ploidy_ht = ploidy_ht.annotate(**sex_ht[ploidy_ht.key])
        ploidy_ht = ploidy_ht.annotate_globals(f_stat_cutoff=f_stat_cutoff)

    if infer_karyotype:
        karyotype_ht = infer_sex_karyotype(
            ploidy_ht, f_stat_cutoff, use_gaussian_mixture_model
        )
        ploidy_ht = ploidy_ht.annotate(**karyotype_ht[ploidy_ht.key])
        ploidy_ht = ploidy_ht.annotate_globals(**karyotype_ht.index_globals())

    return ploidy_ht


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
            ).checkpoint(str(f'{name}_checkpoint.vds'))
            logging.info(f'count post {name} filter:{vds.variant_data.count()}')

    # Infer sex (adds row fields: is_female, var_data_chr20_mean_dp, sex_karyotype)
    sex_ht = annotate_sex(
        vds,
        tmp_prefix=str('annotate_sex'),
        overwrite=not get_config()['workflow'].get('check_intermediates'),
        included_intervals=calling_intervals_ht,
        gt_expr='LGT',
        variants_only_x_ploidy=True,
        variants_only_y_ploidy=False,
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
