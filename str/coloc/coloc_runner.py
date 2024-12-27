#!/usr/bin/env python3

"""
This script performs SNP-only colocalisation analysis betweeen eGenes identified by pseudobulk STR analysis and GWAS signals.
Assumes that the SNP GWAS data has been pre-processed with the following columns: 'chromosome', 'position' (hg38 bp), 'snp'(chromosome_position_refallele_effectallele), 'beta', 'varbeta'

1) Identify eGenes where at least one STR has pval < 5e-8
2) Extract the SNP GWAS data for the cis-window (gene +/- 100kB)
3) Run coloc for each eGene (if the SNP GWAS data has at least one variant with pval <5e-8)
4) Write the results to a TSV file

analysis-runner --dataset "bioheart" \
    --description "Run coloc for eGenes identified by STR analysis" \
    --access-level "test" \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "tenk10k/str/associatr/final_freeze" \
    coloc_runner.py \
    --snp-gwas-file=gs://cpg-bioheart-test/str/Trujillo_methylation_eQTLs/hg38_STRs_SNVs_parsed.tsv \
    --pheno-output-name="Trujillo_methylation_eQTLs" \
    --celltypes='ILC'

"""

import click
import pandas as pd

import hailtop.batch as hb

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, image_path, output_path


def coloc_runner(gwas, probe, eqtl_file_path, celltype, pheno_output_name):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    from cpg_utils.hail_batch import output_path

    ro.r('library(coloc)')
    ro.r('library(tidyverse)')

    with (ro.default_converter + pandas2ri.converter).context():
        gwas_r = ro.conversion.get_conversion().py2rpy(gwas)
    ro.globalenv['gwas_r'] = gwas_r
    ro.r(
        '''
    gwas_r = gwas_r %>% select(beta, varbeta, position,snp)
    gwas_r = gwas_r %>% filter((beta!=0) | (varbeta!=0))
    gwas_r = gwas_r %>% distinct(snp, .keep_all = TRUE)
    gwas_r = gwas_r%>% as.list()
    gwas_r$type = 'quant'
    gwas_r$sdY = 1
    ''',
    )
    eqtl = pd.read_csv(
        eqtl_file_path,
        sep='\t',
    )
    eqtl['beta'] = eqtl['coeff_meta']
    eqtl['se'] = eqtl['se_meta']
    eqtl['position'] = eqtl['pos']
    eqtl['snp'] = eqtl['chr'] + '_' + eqtl['position'].astype(str) + '_' + eqtl['motif']
    eqtl['snp'] = eqtl['snp'].str.replace('-', '_', regex=False)
    gene = eqtl_file_path.split('/')[-1].split('_')[0]
    with (ro.default_converter + pandas2ri.converter).context():
        eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
    ro.globalenv['eqtl_r'] = eqtl_r
    ro.globalenv['gene'] = gene
    ro.r(
        '''
    eqtl_r = eqtl_r %>% filter(!is.na(beta))
    eqtl_r = eqtl_r %>% distinct(snp, .keep_all = TRUE)
    eqtl_r$varbeta = eqtl_r$se**2
    eqtl_r$position = eqtl_r$pos
    eqtl_r = eqtl_r %>% select(beta, varbeta, position, snp)

    eqtl_r = eqtl_r %>% as.list()
    eqtl_r$type = 'quant'
    eqtl_r$sdY = 1


    my.res <- coloc.abf(dataset1=gwas_r,
                    dataset2=eqtl_r)

    p_df <- data.frame(gene,my.res$summary[1], my.res$summary[2], my.res$summary[3], my.res$summary[4], my.res$summary[5], my.res$summary[6])
    names(p_df) <- c('gene', 'nsnps_coloc_tested','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
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
        output_path(f"coloc-tr-snp/sig_str_and_gwas_hit/{pheno_output_name}/{celltype}/{gene}_{probe}_100kb.tsv", 'analysis'),
        sep='\t',
        index=False,
    )


@click.option(
    '--egenes-file',
    help='Path to the eGenes file with FINEMAP and SUSIE probabilities',
    default='gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/finemapped_etrs.csv',
)
@click.option(
    '--snp-cis-dir',
    help='Path to the directory containing the SNP cis results',
    default='gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/snps_and_strs/bioheart_n975_and_tob_n950/meta_results',
)
@click.option(
    '--snp-gwas-file',
    help='Path to the SNP GWAS file',
    default='gs://cpg-bioheart-test/str/gwas_catalog/gcst/gcst-gwas-catalogs/GCST011071_parsed.tsv',
)
@click.option('--celltypes', help='Cell types to run', default='ASDC')
@click.option('--max-parallel-jobs', help='Maximum number of parallel jobs to run', default=500)
@click.option('--pheno-output-name', help='Phenotype output name', default='covid_GCST011071')
@click.option('--job-cpu', help='Number of CPUs to use for each job', default=0.25)
@click.command()
def main(snp_cis_dir, egenes_file, celltypes, snp_gwas_file, pheno_output_name, max_parallel_jobs, job_cpu):
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
    var_table = pd.read_csv(
        'gs://cpg-bioheart-test/tenk10k/saige-qtl/300libraries_n1925_adata_raw_var.csv',
    )
    hg38_map = pd.read_csv(
        snp_gwas_file,
        sep='\t',
    )

    hg38_map['p_value'] = hg38_map['p_value'].astype(float)
    hg38_map['position'] = hg38_map['position'].astype(int)

    # read in eGenes file
    egenes = pd.read_csv(
        egenes_file)

    result_df_cfm_str = egenes
    #result_df_cfm['variant_type'] = result_df_cfm['motif'].str.contains('-').map({True: 'SNV', False: 'STR'})
    #result_df_cfm_str = result_df_cfm[result_df_cfm['variant_type'] == 'STR']  # filter for STRs
    result_df_cfm_str = result_df_cfm_str[result_df_cfm_str['pval_meta'] < 5e-8]  # filter for STRs with p-value < 5e-8
    result_df_cfm_str = result_df_cfm_str.drop_duplicates(
        subset=['gene_name', 'cell_type'],
    )  # drop duplicates (ie pull out the distinct genes in each celltype)
    result_df_cfm_str['gene'] = result_df_cfm_str['gene_name']
    result_df_cfm_str['celltype'] = result_df_cfm_str['cell_type']
    b = get_batch(name=f'Run coloc:{pheno_output_name}')

    for celltype in celltypes.split(','):
        result_df_cfm_str_celltype = result_df_cfm_str[
            result_df_cfm_str['celltype'] == celltype
        ]  # filter for the celltype of interest
        for gene in result_df_cfm_str_celltype['gene']:
            chrom = result_df_cfm_str_celltype[result_df_cfm_str_celltype['gene'] == gene]['chr'].iloc[0]

            if to_path(f'{snp_cis_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv').exists():
                print('Cis results for ' + gene + ' exist: proceed with coloc')

                # extract the coordinates for the cis-window (gene +/- 100kB)
                gene_table = var_table[var_table['gene_ids'] == gene]
                start = float(gene_table['start'].astype(float)) - 100000
                end = float(gene_table['end'].astype(float)) + 100000
                chrom = gene_table['chr'].iloc[0]
                hg38_map_chr = hg38_map[hg38_map['chromosome'] == (chrom)]
                hg38_map_chr_start = hg38_map_chr[hg38_map_chr['position'] >= start]
                hg38_map_chr_start_end = hg38_map_chr_start[hg38_map_chr_start['position'] <= end]
                for probe in hg38_map_chr_start_end['ProbeID'].unique():
                    if to_path(
                        output_path(
                            f"coloc-tr-snp/sig_str_and_gwas_hit/{pheno_output_name}/{celltype}/{gene}_{probe}_100kb.tsv",
                            'analysis',
                        ),
                    ).exists():
                        continue
                    hg38_map_chr_start_end_probe = hg38_map_chr_start_end[hg38_map_chr_start_end['ProbeID'] == probe]
                    if hg38_map_chr_start_end_probe.empty:
                        print('No SNP GWAS data for ' + gene + f' {probe} combination' +' in the cis-window: skipping....')
                        continue
                    # check if the p-value column contains at least one value which is <5e-8:
                    if hg38_map_chr_start_end_probe['p_value'].min() >= 5e-8:
                        print('No significant SNP GWAS data for ' + gene + f' {probe} combination' + ' in the cis-window: skipping....')
                        continue
                    print('Extracted SNP GWAS data for ' +gene + f' {probe} combination')

                    # run coloc
                    coloc_job = b.new_python_job(
                        f'Coloc for {gene} and {probe}: {celltype}',
                    )
                    coloc_job.image(image_path('r-meta'))
                    coloc_job.cpu(job_cpu)
                    coloc_job.call(
                        coloc_runner,
                        hg38_map_chr_start_end_probe,
                        probe,
                        f'{snp_cis_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv',
                        celltype,
                        pheno_output_name,
                    )
                    manage_concurrency_for_job(coloc_job)

            else:
                print('No cis results for ' + gene + ' exist: skipping....')

    b.run(wait=False)


if __name__ == '__main__':
    main()
