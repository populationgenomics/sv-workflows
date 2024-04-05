#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda
"""
This script applies filters to a ExpansionHunter MT, and outputs chromosome-specific VCF ready for input into associaTR.
Input MT should be the output of qc_annotator.py

Applied filters:
- remove monomorphic variants
- locus-level callrate >=0.9
- sample-level call rate >=0.99
- observed heterozygosity >=0.00995 (corresponds to a MAF = 0.5%)
- locus-level HWEP (binom definition) >=10^-6
- set calls outside of [-30,20] relative to mode to NA
- enforce locus-level call rate >=0.9


 analysis-runner --dataset "bioheart" \
    --description "Hail QC for associaTR" \
    --access-level "test" \
    --output-dir "str/associatr/input_files" \
    qc_filters_cohort_spec.py --mt-path=gs://cpg-bioheart-test/str/polymorphic_run_n2045/annotated_mt/v2/str_annotated.mt

"""

import hail as hl
import click
from ast import literal_eval
from cpg_utils import to_path


from cpg_utils.config import get_config

from cpg_utils.hail_batch import get_batch, init_batch, output_path

config = get_config()

BCFTOOLS_IMAGE = config['images']['bcftools']


def qc_filter(mt_path, version):
    """
    Applies QC filters to the input MT
    """
    init_batch(worker_memory='highmem')

    # read in mt
    mt = hl.read_matrix_table(mt_path)

    # remove monomorphic variants, set locus level call rate >=0.9, observed heterozygosity >=0.00995, locus level HWEP (binom definition) >=10^-6
    mt = mt.filter_rows(
        (mt.num_alleles > 1)
        & (mt.variant_qc.call_rate >= 0.9)
        & (mt.obs_het >= 0.00995)
        & (mt.binom_hwep >= 0.000001)
    )

    # set sample level call rate >=0.99
    mt = mt.filter_cols(mt.sample_qc.call_rate >= 0.99)

    # set big expansions/deletions beyond [-30,20] relative to mode allele to NA
    condition_allele_1 = (mt.allele_1_minus_mode > 20) | (mt.allele_1_minus_mode < -30)
    mt = mt.annotate_entries(
        GT=hl.if_else(
            condition_allele_1,
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            condition_allele_1,
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            condition_allele_1,
            hl.missing('str'),
            mt.SO,
        ),
        allele_1_rep_length=hl.if_else(
            condition_allele_1,
            hl.missing('int32'),
            mt.allele_1_rep_length,
        ),
    )

    condition_allele_2 = (mt.allele_2_minus_mode > 20) | (mt.allele_2_minus_mode < -30)

    mt = mt.annotate_entries(
        GT=hl.if_else(
            condition_allele_2,
            hl.missing('call'),
            mt.GT,
        ),
        ADFL=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADFL,
        ),
        ADIR=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADIR,
        ),
        ADSP=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.ADSP,
        ),
        LC=hl.if_else(
            condition_allele_2,
            hl.missing('float64'),
            mt.LC,
        ),
        REPCI=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.REPCI,
        ),
        REPCN=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.REPCN,
        ),
        SO=hl.if_else(
            condition_allele_2,
            hl.missing('str'),
            mt.SO,
        ),
        allele_2_rep_length=hl.if_else(
            condition_allele_2,
            hl.missing('int32'),
            mt.allele_2_rep_length,
        ),
    )

    # calculate proportion of GTs that are defined per locus (after applying call-level filters, variant_qc.call_rate is not accurate anymore)
    mt = mt.annotate_rows(
        prop_GT_exists=hl.agg.count_where(hl.is_defined(mt.GT))
        / (mt.variant_qc.n_called + mt.variant_qc.n_not_called)
    )
    # re-enforce locus level call rate >=0.9
    mt = mt.filter_rows(mt.prop_GT_exists >= 0.9)

    ### APPLY SUBSETS FOR COHORT SPECIFIC MTS
    # remove chrX from analysis
    mt = mt.filter_rows((hl.str(mt.locus.contig).startswith('chrX')), keep=False)
    table = hl.import_table(
        'gs://cpg-bioheart-test/str/sample-sex-mapping/sample_karyotype_sex_mapping.csv',
        delimiter=',',
        impute=True,
    )

    # Define a function to determine the cohort based on the 'id' column
    def get_cohort(id):
        return hl.if_else(
            hl.str(id).startswith('CT'),
            'bioheart',
            hl.if_else(hl.str(id).startswith('TOB'), 'tob', 'Unknown'),
        )

    # Add a new column 'cohort' based on the 'id' column using the defined function
    table = table.annotate(cohort=get_cohort(table.external_id))
    table = table.key_by('s')
    mt = mt.annotate_cols(cohort=table[mt.s].cohort)

    table_geno_pcs = hl.import_table(
        'gs://cpg-bioheart-test/str/anndata/saige-qtl/input_files/covariates/sex_age_geno_pcs_tob_bioheart.csv',
        delimiter=',',
        impute=True,
    )

    table_geno_pcs = table_geno_pcs.key_by('sample_id')
    mt = mt.annotate_cols(geno_pc1=hl.float(table_geno_pcs[mt.s].geno_PC1))
    mt = mt.annotate_cols(geno_pc6=hl.float(table_geno_pcs[mt.s].geno_PC6))
    # remove ancestry outliers
    mt = mt.filter_cols(
        (mt.geno_pc1 >= -0.05) & (mt.geno_pc6 <= 0.05) & (mt.geno_pc6 >= -0.05)
    )
    with to_path(
        'gs://cpg-bioheart-test/str/associatr/input_files/remove-samples.txt'
    ).open() as f:
        array_string = f.read().strip()
        remove_samples = literal_eval(array_string)

    # remove related individuals
    mt = mt.filter_cols(hl.literal(remove_samples).contains(mt.s), keep=False)

    ### calculate col specific values and export
    mt = mt.annotate_cols(
        col_sum_allele_1_is_not_mode=hl.agg.sum(hl.cond(mt.allele_1_is_non_mode, 1, 0)),
        col_sum_allele_2_is_not_mode=hl.agg.sum(hl.cond(mt.allele_2_is_non_mode, 1, 0)),
    )
    mt = mt.annotate_cols(
        col_sum_alleles_is_not_mode=mt.col_sum_allele_1_is_not_mode
        + mt.col_sum_allele_2_is_not_mode
    )
    mt = mt.annotate_cols(missing_calls_count=hl.agg.count_where(hl.is_missing(mt.GT)))

    mt.cols().write('gs://cpg-bioheart-test/str/batch_debug/cols.ht', overwrite=True)

    for cohort in ['tob', 'bioheart']:
        mt_cohort = mt.filter_cols(mt.cohort == cohort)

        # annotate rows with counts of alleles that appear in each VARID
        ht = mt_cohort.select_rows(
            alleles_rep_lengths=hl.agg.collect(mt_cohort.allele_1_rep_length).extend(
                hl.agg.collect(mt_cohort.allele_2_rep_length)
            )
        ).rows()

        # explode the allele_array to create one row for each element in the array
        exploded_table = ht.explode(ht.alleles_rep_lengths)

        # aggregate the counts for each distinct value in the exploded allele_array
        aggregated_table = exploded_table.group_by(exploded_table.REPID).aggregate(
            allele_array_counts=hl.agg.counter(exploded_table.alleles_rep_lengths)
        )

        # annotate the mt with the aggregated counts
        mt_cohort = mt_cohort.annotate_rows(
            aggregated_info_cohort=aggregated_table[mt_cohort.REPID]
        )
        mt_cohort = mt_cohort.annotate_rows(
            cohort_num_alleles=hl.len(
                mt_cohort.aggregated_info_cohort.allele_array_counts.values()
            )
        )
        mt_cohort = mt_cohort.annotate_rows(
            sum_allele_1_is_not_mode=hl.agg.sum(
                hl.cond(mt_cohort.allele_1_is_non_mode, 1, 0)
            ),
            sum_allele_2_is_not_mode=hl.agg.sum(
                hl.cond(mt_cohort.allele_2_is_non_mode, 1, 0)
            ),
        )
        mt_cohort = mt_cohort.annotate_rows(
            sum_alleles_is_not_mode=mt_cohort.sum_allele_1_is_not_mode
            + mt_cohort.sum_allele_2_is_not_mode
        )

        # proportion of alleles that are not the mode allele
        mt_cohort = mt_cohort.annotate_rows(
            prop_alleles_is_not_mode=mt_cohort.sum_alleles_is_not_mode
            / hl.sum(mt_cohort.aggregated_info_cohort.allele_array_counts.values())
        )

        mt_cohort.rows().write(
            f'gs://cpg-bioheart-test/str/batch_debug/rows_{cohort}.ht', overwrite=True
        )


@click.option(
    '--mt-path',
    help='GCS file path to mt (output of qc_annotator.py)',
    type=str,
)
@click.option('--hail-storage', help='Hail storage', type=str, default='0G')
@click.option('--hail-cpu', help='Hail CPU', type=int, default=1)
@click.option('--hail-memory', help='Hail memory', type=str, default='standard')
@click.option('--version', help='version of the output files', type=str, default='v1')
@click.command()
def main(
    mt_path,
    hail_storage,
    hail_cpu,
    hail_memory,
    version,
):
    """
    Runner to apply QC filters to input MT, and bgzip and tabix.
    """

    b = get_batch('Apply QC filters to MT and export chr-specific VCFs')

    hail_job = b.new_python_job(name=f'QC filters')
    hail_job.image(config['workflow']['driver_image'])
    hail_job.storage(hail_storage)
    hail_job.cpu(hail_cpu)
    hail_job.memory(hail_memory)
    hail_job.call(qc_filter, mt_path, version)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
