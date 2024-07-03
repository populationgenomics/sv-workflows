#!/usr/bin/env python3

"""
This script performs colocalisation analysis betweeen eGenes identified by pseudobulk STR analysis and GWAS signals from UKBB Gymrek STR catalogs.
1) Identify eGenes (FDR <5%) from the STR analysis
2) Extract the STR GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene (if STR eQTL data is available)
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by STR analysis" \
    --access-level "test" \
    --output-dir "str/associatr" \
    coloc_gymrek_ukbb_runner.py \
    --celltypes "ASDC" \
    --egenes-dir='gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/all_cell_types_all_genes_sig_only.tsv' \
    --pheno "albumin"
"""

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils.hail_batch import get_batch, image_path


def cyclical_shifts(s):
    return [s[i:] + s[:i] for i in range(len(s))]


def coloc_runner(result_df_cfm_str_celltype, var_table, str_gwas_file, snp_gwas_file, eqtl_cis_dir, celltype, pheno):
    import gzip

    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path

    with gzip.open(to_path(str_gwas_file), 'rb') as f:
            str_gwas = pd.read_csv(
                f,
                sep='\t',
                usecols=[
                    'chromosome',
                    'beta',
                    'standard_error',
                    'p_value',
                    'repeat_unit',
                    'start_pos (hg38)',
                    'end_pos (hg38)',
                ],
            )
    with gzip.open(to_path(snp_gwas_file), 'rb') as snp_f:
        snp_gwas = pd.read_csv(
            snp_f,
            sep='\t',
        )

    for gene in result_df_cfm_str_celltype['gene']:
        if to_path(
            output_path(
                f"coloc/sig_str_and_gwas_hit/gymrek-ukbb-{pheno}/{celltype}/{gene}_100kb.tsv",
                'analysis',
            ),
        ).exists():
            print('Coloc results already processed for ' + gene + ': skipping....')
            continue
        # extract the coordinates for the cis-window (gene +/- 100kB)
        gene_table = var_table[var_table['gene_ids'] == gene]
        start = float(gene_table['start'].astype(float)) - 100000
        end = float(gene_table['end'].astype(float)) + 100000
        chrom = gene_table['chr'].iloc[0][3:]

        ## subset the STR GWAS data for the cis-window
        str_gwas_subset = str_gwas[
            (str_gwas['chromosome'] == int(chrom))
            & (str_gwas['start_pos (hg38)'] >= start)
            & (str_gwas['start_pos (hg38)'] <= end)
        ]

        ## subset the SNP GWAS data for the cis-window
        snp_gwas_subset = snp_gwas[
            (snp_gwas['chromosome'] == gene_table['chr'].iloc[0])
            & (snp_gwas['position'] >= start)
            & (snp_gwas['position'] <= end)
        ]

        if str_gwas_subset.empty and snp_gwas_subset.empty:
            print('No GWAS data for ' + gene + ' in the cis-window: skipping....')
            continue

        if str_gwas_subset['p_value'].min() > 5e-8 and snp_gwas_subset['p_value'].min() > 5e-8:
            print('No significant GWAS data for ' + gene + ' in the cis-window: skipping....')
            continue

        print('Extracted GWAS data for ' + gene)  # proceed to coloc for this gene

        eqtl_file_path = f'{eqtl_cis_dir}/{celltype}/chr{chrom}/{gene}_100000bp_meta_results.tsv'
        eqtl = pd.read_csv(eqtl_file_path, sep='\t')

        # create an empty df with the following columns: 'chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value'
        gwas_str_harmonised = pd.DataFrame(columns=['chromosome', 'position', 'varbeta', 'beta', 'snp', 'p_value'])

        # harmonise Gymrek STR definitons with ours
        eqtl_str = eqtl[~eqtl['motif'].str.contains('-')]  # str entries have no '-' in motif column
        for index, gwas_row in str_gwas_subset.iterrows():
            for index2, eqtl_row in eqtl_str.iterrows():
                if (
                    eqtl_row['pos'] >= (gwas_row['start_pos (hg38)'] - 1)
                    and eqtl_row['pos'] <= gwas_row['end_pos (hg38)']
                ):
                    if eqtl_row['motif'] in cyclical_shifts(gwas_row['repeat_unit']):
                        new_entry = pd.DataFrame(
                            [
                                {
                                    'chromosome': eqtl_row["chr"],
                                    'position': eqtl_row['pos'],
                                    'varbeta': gwas_row['standard_error'] ** 2,
                                    'beta': gwas_row['beta'],
                                    'snp': f'{eqtl_row["chr"]}_{eqtl_row["pos"]}_{eqtl_row["motif"]}',
                                    'p_value': gwas_row['p_value'],
                                },
                            ],
                        )
                        gwas_str_harmonised = pd.concat([gwas_str_harmonised, new_entry], ignore_index=True)

                        continue

        # concatenate gwas_str with gwas_snp (row wise)
        gwas = pd.concat([gwas_str_harmonised, snp_gwas_subset], ignore_index=True)
        gwas = gwas.sort_values(by=['snp', 'p_value'])
        # Drop duplicate variants (keep the one with the lowest p-value)
        gwas = gwas.drop_duplicates(subset='snp', keep='first')

        gwas = gwas[~gwas['beta'].isna()]

        ro.r('library(coloc)')
        ro.r('library(tidyverse)')

        with (ro.default_converter + pandas2ri.converter).context():
            gwas_r = ro.conversion.get_conversion().py2rpy(gwas)
        print('loaded in gwas_r')
        ro.globalenv['gwas_r'] = gwas_r
        ro.r(
            '''
                gwas_r = gwas_r %>% distinct(beta, varbeta, position,snp)
                gwas_r = gwas_r%>% as.list()
                gwas_r$type = 'quant'
                gwas_r$sdY = 1

                ''',
        )
        eqtl['beta'] = eqtl['coeff_meta']
        eqtl['varbeta'] = eqtl['se_meta'] ** 2
        eqtl['position'] = eqtl['pos']
        eqtl['snp'] = eqtl['chr'] + '_' + eqtl['position'].astype(str) + '_' + eqtl['motif']
        eqtl['snp'] = eqtl['snp'].str.replace('-', '_', regex=False)

        eqtl = eqtl.sort_values(by=['snp', 'pval_meta'])
        # Drop duplicate variants (keep the one with the lowest p-value)
        eqtl = eqtl.drop_duplicates(subset='snp', keep='first')

        with (ro.default_converter + pandas2ri.converter).context():
            eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
        print('loaded in eqtl_r')
        ro.globalenv['eqtl_r'] = eqtl_r
        ro.globalenv['gene'] = gene
        ro.r(
            '''
                eqtl_r = eqtl_r%>% distinct(beta, varbeta, position, snp)
                eqtl_r = eqtl_r %>% as.list()
                eqtl_r$type = 'quant'
                eqtl_r$sdY = 1

                p4 <- tryCatch({

                my.res <- coloc.abf(dataset1=gwas_r,
                                dataset2=eqtl_r)

                p_df <- data.frame(gene,my.res$summary[1], my.res$summary[2], my.res$summary[3], my.res$summary[4], my.res$summary[5], my.res$summary[6])
                names(p_df) <- c('gene', 'nvariants_coloc_tested','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')


                }, error = function(e) {

                print(paste("An error occurred:", e))
                0
                })

                ''',
        )

        # convert to pandas df
        with (ro.default_converter + pandas2ri.converter).context():
            pd_p4_df = ro.conversion.get_conversion().rpy2py(ro.r('p_df'))

        # add cell type and chrom annotation to df
        pd_p4_df['celltype'] = celltype
        pd_p4_df['chrom'] = eqtl['chr'].iloc[0]

        # write to GCS
        pd_p4_df.to_csv(
            output_path(f"coloc/sig_str_and_gwas_hit/gymrek-ukbb-{pheno}/{celltype}/{gene}_100kb.tsv", 'analysis'),
            sep='\t',
            index=False,
        )


@click.option(
    '--egenes-dir',
    help='Path to the eGenes dir',
    default='gs://cpg-bioheart-test-analysis/str/associatr/tob_n1055_and_bioheart_n990/DL_random_model/meta_results/fdr_qvals/using_acat',
)
@click.option(
    '--eqtl-cis-dir',
    help='Path to the directory containing the STR+SNP eQTL cis results',
    default='gs://cpg-bioheart-test/str/associatr/snps_and_strs/tob_n1055_and_bioheart_n990/meta_results/meta_results',
)
@click.option('--celltypes', help='Cell type for which the eGenes were identified', default='CD4_TCM')
@click.option(
    '--var-annotation-file',
    help='Gene annotation file path',
    default='gs://cpg-bioheart-test/str/240_libraries_tenk10kp1_v2/concatenated_gene_info_donor_info_var.csv',
)
@click.option('--pheno', help='Phenotype to use for coloc', default='alanine_aminotransferase')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs', default=500)
@click.command()
def main(eqtl_cis_dir, egenes_dir, celltypes, var_annotation_file, pheno, max_parallel_jobs):
    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= max_parallel_jobs:
            job.depends_on(_dependent_jobs[-max_parallel_jobs])
        _dependent_jobs.append(job)

    # read in gene annotation file
    var_table = pd.read_csv(var_annotation_file)

    for phenotype in pheno.split(','):
        str_gwas_file = f'gs://cpg-bioheart-test/str/gymrek-ukbb-str-gwas-catalogs/gymrek-ukbb-str-gwas-catalogs/white_british_{phenotype}_str_gwas_results.tab.gz'
        snp_gwas_file = f'gs://cpg-bioheart-test/str/gymrek-ukbb-snp-gwas-catalogs/white_british_{phenotype}_snp_gwas_results_hg38.tab.gz'

        # read in eGenes file
        egenes = pd.read_csv(
            egenes_dir,
            sep='\t',
            usecols=['chr', 'pos', 'pval_meta', 'motif', 'susie_pip', 'gene', 'finemap_prob', 'celltype', 'ref_len'],
        )

        result_df_cfm = egenes
        result_df_cfm['variant_type'] = result_df_cfm['motif'].str.contains('-').map({True: 'SNV', False: 'STR'})
        result_df_cfm_str = result_df_cfm[result_df_cfm['variant_type'] == 'STR']  # filter for STRs
        result_df_cfm_str = result_df_cfm_str[
            result_df_cfm_str['pval_meta'] < 5e-8
        ]  # filter for STRs with p-value < 5e-8
        result_df_cfm_str = result_df_cfm_str.drop_duplicates(
            subset=['gene', 'celltype'],
        )  # drop duplicates (ie pull out the distinct genes in each celltype)
        result_df_cfm_str['gene'] = result_df_cfm_str['gene'].str.replace(
            '.tsv',
            '',
            regex=False,
        )  # remove .tsv from gene names (artefact of the data file)

        b = get_batch(name='Run coloc')

        for celltype in celltypes.split(','):
            result_df_cfm_str_celltype = result_df_cfm_str[
                result_df_cfm_str['celltype'] == celltype
            ]  # filter for the celltype of interest

            # run one coloc job for every celltype-phenotype combination
            coloc_job = b.new_python_job(
                f'Coloc for {celltype} and {phenotype}',
            )
            coloc_job.cpu(8)
            coloc_job.image(image_path('r-meta'))
            coloc_job.call(
                coloc_runner,
                result_df_cfm_str_celltype,
                var_table,
                str_gwas_file,
                snp_gwas_file,
                eqtl_cis_dir,
                celltype,
                phenotype,
            )
            manage_concurrency_for_job(coloc_job)

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()
