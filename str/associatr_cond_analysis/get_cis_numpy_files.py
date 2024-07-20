#!/usr/bin/env python3
# pylint: disable=too-many-arguments,too-many-locals

"""
This script aims to prepare inputs for conditional analysis by:
 - loading in the list of eGenes and filtering for loci where the STR is the lead signal and subsequently:
 - output gene lists for each cell type and chromosome (after filtering out lowly expressed genes)
 - output cis window files for each gene with scRNA data (cell type + chr specific)
 - optionally removes samples based on a provided sample file
 - perform rank-based inverse normal transformation on pseudobulk data (per gene basis)
 - output gene-level phenotype and covariate numpy objects for input into associatr, with lead STR genotypes as a covariate.

 analysis-runner  --config get_cis_numpy_files.toml --dataset "bioheart" --access-level "test" \
--description "get cis and numpy" --output-dir "str/associatr/cond_analysis/tob_n1055" \
--image australia-southeast1-docker.pkg.dev/cpg-common/images/scanpy:1.9.3 \
python3 get_cis_numpy_files.py

"""
import json
from ast import literal_eval

import pandas as pd
from cyvcf2 import VCF

import hail as hl
import hailtop.batch as hb

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, image_path, init_batch


def extract_snp_genotypes(vcf_file, loci):
    """
    Helper function to extract genotypes (SNPs) from a VCF file; target loci specified as a list (can be single or multiple)

    """
    # Read the VCF file
    vcf_reader = VCF(vcf_file)

    results = pd.DataFrame()
    results['sample_id'] = vcf_reader.samples

    # Iterate through the records in the VCF file for each locus
    for locus in loci:
        chrom, pos = locus.split(':')
        pos = int(pos)

        for record in vcf_reader(f'{chrom}:{pos}-{pos}'):
            if record.CHROM == chrom and record.POS == pos:
                gt = record.gt_types
                gt[gt == 3] = 2  # HOM ALT is coded as 3; change it to 2
                results[locus] = gt
                break

    # Convert results to a DataFrame
    return results


def extract_str_genotypes(vcf_file, loci, motifs):
    """
    Helper function to extract genotypes (STRs) from a VCF file; target loci specified as a list (can be single or multiple)

    """
    import numpy as np
    import pandas as pd
    from cyvcf2 import VCF

    # Read the VCF file
    vcf_reader = VCF(vcf_file)

    results = pd.DataFrame()
    results['sample_id'] = vcf_reader.samples

    for locus, motif in zip(loci, motifs):
        chrom, pos = locus.split(':')
        pos = int(pos)

        for record in vcf_reader(f'{chrom}:{pos}-{pos}'):
            # check motif as well because two STRs can have same coord but different motif
            if record.CHROM == chrom and record.POS == pos and record.INFO.get('RU') == motif:
                genotypes = record.format('REPCN')
                # Replace '.' with '-99/-99' to handle missing values
                genotypes = np.where(genotypes == '.', '-99/-99', genotypes)

                # Split each element by '/'
                split_genotypes = [genotype.split('/') for genotype in genotypes]

                # Convert split_genotypes into a numpy array for easier manipulation
                split_genotypes_array = np.array(split_genotypes)

                # Convert the strings to integers and sum them row-wise
                sums = np.sum(split_genotypes_array.astype(int), axis=1)
                # set dummy -198 value to np.nan
                sums = np.where(sums == -198, np.nan, sums)

                results[locus] = sums
                break

    return results


def cis_window_numpy_extractor(
    egenes_path,
    input_h5ad_dir,
    input_pseudobulk_dir,
    input_cov_dir,
    chromosome,
    cell_type,
    cis_window,
    version,
    chrom_len,
    min_pct,
    remove_samples_file,
    str_input,
):
    """
    Creates gene-specific cis window files and phenotype-covariate numpy objects

    """
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy.stats import norm

    import hail as hl

    from cpg_utils import to_path
    from cpg_utils.config import get_config
    from cpg_utils.hail_batch import init_batch, output_path

    init_batch()

    chromosome = f'chr{chromosome}'
    egenes = pd.read_csv(egenes_path)
    egenes_cell = egenes[egenes['cell_type'] == cell_type]

    egenes_cell_chrom = egenes_cell[egenes_cell['chr'] == chromosome]

    # read in anndata object because anndata.vars has the start, end coordinates of each gene
    h5ad_file_path = f'{input_h5ad_dir}/{cell_type}_{chromosome}.h5ad'
    expression_h5ad_path = to_path(h5ad_file_path).copy('here.h5ad')
    adata = sc.read_h5ad(expression_h5ad_path)

    # read in pseudobulk and covariate files
    pseudobulk_path = f'{input_pseudobulk_dir}/{cell_type}/{cell_type}_{chromosome}_pseudobulk.csv'
    pseudobulk = pd.read_csv(pseudobulk_path)
    covariate_path = f'{input_cov_dir}/{cell_type}_covariates.csv'
    covariates = pd.read_csv(covariate_path)

    # extract genes in egenes df
    gene_names = egenes_cell_chrom['gene_name']

    for gene in gene_names:
        try:
            eqtl_results = pd.read_csv(
                f'{get_config()["get_cis_numpy"]["snp_str_meta"]}/{cell_type}/{chromosome}/{gene}_100000bp_meta_results.tsv',
                sep='\t',
            )
        except FileNotFoundError:
            print(f'No eQTL results found for {gene}... skipping')
            continue
        # get row(s) with minimum p-value
        min_pval = eqtl_results['pval_meta'].min()
        smallest_pval_rows = eqtl_results[eqtl_results['pval_meta'] == min_pval]
        # check if all rows are SNPs:
        all_motif_dash = smallest_pval_rows['motif'].str.contains('-').all()
        if all_motif_dash:
            print(f'Lead signal(s) is a SNP for {gene}... skipping')
            continue

        # get characteristics of the lead STR
        lead_str = smallest_pval_rows[smallest_pval_rows['motif'] != '-']
        lead_str_motif = [lead_str['motif'].values[0]]
        lead_str_locus = [f'{chromosome}:{lead_str["pos"].values[0]}']

        # get gene body position (start and end) and add window
        start_coord = adata.var[adata.var.index == gene]['start']
        end_coord = adata.var[adata.var.index == gene]['end']

        left_boundary = max(1, int(start_coord.iloc[0]) - int(cis_window))
        right_boundary = min(int(end_coord.iloc[0]) + int(cis_window), chrom_len)

        data = {'chromosome': chromosome, 'start': left_boundary, 'end': right_boundary}
        ofile_path = output_path(
            f'input_files/cis_window_files/{version}/{cell_type}/{chromosome}/{gene}_{cis_window}bp.bed',
        )
        # write cis window file to gcp
        pd.DataFrame(data, index=[gene]).to_csv(ofile_path, sep='\t', header=False, index=False)

        # make the phenotype-covariate numpy objects
        pseudobulk.rename(columns={'individual': 'sample_id'}, inplace=True)  # noqa: PD002
        gene_pheno = pseudobulk[['sample_id', gene]]

        # remove samples that are in the remove_samples_file
        if remove_samples_file:
            with to_path(remove_samples_file).open() as f:
                array_string = f.read().strip()
                remove_samples = literal_eval(array_string)
                gene_pheno = gene_pheno[~gene_pheno['sample_id'].isin(remove_samples)]

        # rank-based inverse normal transformation based on R's orderNorm()
        # Rank the values
        gene_pheno.loc[:, 'gene_rank'] = gene_pheno[gene].rank()
        # Calculate the percentile of each rank
        gene_pheno.loc[:, 'gene_percentile'] = (gene_pheno.loc[:, 'gene_rank'] - 0.5) / (len(gene_pheno))
        # Use the inverse normal cumulative distribution function (quantile function) to transform percentiles to normal distribution values
        gene_pheno.loc[:, 'gene_inverse_normal'] = norm.ppf(gene_pheno.loc[:, 'gene_percentile'])
        gene_pheno = gene_pheno[['sample_id', 'gene_inverse_normal']]

        gene_pheno_cov = gene_pheno.merge(covariates, on='sample_id', how='inner')

        # add STR genotypes we would like to condition on
        if str_input is not None:
            str_genotype_df = extract_str_genotypes(str_input['vcf'], lead_str_locus, lead_str_motif)
            # append 'CPG' prefix to 'id' column
            str_genotype_df['sample_id'] = str_genotype_df['sample_id'].apply(lambda x: 'CPG' + x)
            gene_pheno_cov = gene_pheno_cov.merge(str_genotype_df, on='sample_id', how='inner')

        # filter for samples that were assigned a CPG ID; unassigned samples after demultiplexing will not have a CPG ID
        gene_pheno_cov = gene_pheno_cov[gene_pheno_cov['sample_id'].str.startswith('CPG')]

        gene_pheno_cov['sample_id'] = gene_pheno_cov['sample_id'].str[
            3:
        ]  # remove CPG prefix because associatr expects id to be numeric

        gene_pheno_cov['sample_id'] = gene_pheno_cov['sample_id'].astype(float)

        gene_pheno_cov = gene_pheno_cov.to_numpy()
        with hl.hadoop_open(
            output_path(
                f'input_files/pheno_cov_numpy/{version}/{cell_type}/{chromosome}/{gene}_pheno_cov.npy',
            ),
            'wb',
        ) as f:
            np.save(f, gene_pheno_cov)

        # write filtered gene name to a json file
        with to_path(
            output_path(
                f'input_files/scRNA_gene_lists/{min_pct}_min_pct_cells_expressed/{cell_type}/{chromosome}_{cell_type}_{gene}.json',
            ),
        ).open('w') as write_handle:
            json.dump(gene, write_handle)


def main():
    """
    Run cis window extraction and phenotype/covariate numpy object creation
    """
    b = get_batch(name='get cis_numpy files')

    # Setup MAX concurrency by genes
    _dependent_jobs: list[hb.batch.job.Job] = []

    def manage_concurrency_for_job(job: hb.batch.job.Job):
        """
        To avoid having too many jobs running at once, we have to limit concurrency.
        """
        if len(_dependent_jobs) >= get_config()['get_cis_numpy']['max_parallel_jobs']:
            job.depends_on(_dependent_jobs[-get_config()['get_cis_numpy']['max_parallel_jobs']])
        _dependent_jobs.append(job)

    egenes = pd.read_csv(get_config()['get_cis_numpy']['egene'])
    celltypes = egenes['cell_type'].unique()

    for cell_type in celltypes:
        egenes_cell = egenes[egenes['cell_type'] == cell_type]
        for chrom in range(1, 23):
            egenes_cell_chrom = egenes_cell[egenes_cell['chr'] == f'chr{chrom}']
            if egenes_cell_chrom.empty:
                continue
            init_batch()
            chrom_len = hl.get_reference('GRCh38').lengths[f'chr{chrom}']
            j = b.new_python_job(
                name=f'Extract cis window & phenotype and covariate numpy object for {cell_type}: {chrom}',
            )
            j.image(image_path('scanpy'))
            j.cpu(get_config()['get_cis_numpy']['job_cpu'])
            j.memory(get_config()['get_cis_numpy']['job_memory'])
            j.storage(get_config()['get_cis_numpy']['job_storage'])

            if get_config()['get_cis_numpy']['str_vcf_dir'] != 'gs://cpg-bioheart-test':
                str_vcf_dir = get_config()['get_cis_numpy']['str_vcf_dir']
                str_vcf_path = f'{str_vcf_dir}/hail_filtered_chr{chrom}.vcf.bgz'
                str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'})
            else:  # no STR VCF provided
                str_input = None
            j.call(
                cis_window_numpy_extractor,
                get_config()['get_cis_numpy']['egene'],
                get_config()['get_cis_numpy']['input_h5ad_dir'],
                get_config()['get_cis_numpy']['input_pseudobulk_dir'],
                get_config()['get_cis_numpy']['input_cov_dir'],
                chrom,
                cell_type,
                get_config()['get_cis_numpy']['cis_window'],
                get_config()['get_cis_numpy']['version'],
                chrom_len,
                get_config()['get_cis_numpy']['min_pct'],
                get_config()['get_cis_numpy']['remove_samples_file'],
                str_input,
            )

            manage_concurrency_for_job(j)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,too-many-arguments
