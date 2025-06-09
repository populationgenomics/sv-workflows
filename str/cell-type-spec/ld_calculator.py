#!/usr/bin/env python3

"""
This script calculates linkage disequilibrium (LD) of distinct variants found in the same eGene.
Used for differentiating scenario 3,4,5 in the Cuomo et al. cell typ specificity definitions.

analysis-runner --dataset "tenk10k" --access-level "test"  --description "Prepare LD files"  --output-dir "str/cell-spec-cuomo/v2" ld_calculator.py

"""

import pandas as pd
import pandas as pd
from cpg_utils.hail_batch import get_batch
import click


def build_tr_genotype_matrix(vcf_reader, chrom, gene_start, gene_end, estr_egene, samples_list):
    """

    Builds a genotype matrix for a given egene region from a VCF reader.
    """
    samples = vcf_reader.samples

    dosage_dict = {sample: {} for sample in samples}

    region = f"{chrom}:{gene_start}-{gene_end}"
    print(region)

    for variant in vcf_reader(region):
        variantid = f"{variant.POS}{variant.INFO.get('RU')}"

        if variantid not in set(estr_egene['variantid']):
            continue

        try:
            repcn = variant.format('REPCN')
        except KeyError:
            continue

        for i, sample in enumerate(samples):
            try:
                allele_1, allele_2 = map(int, repcn[i].split('/'))
                dosage = allele_1 + allele_2
            except:
                dosage = None
            dosage_dict[sample][variantid] = dosage

    # Convert to DataFrame: rows = individuals, columns = variant IDs
    genotype_df = pd.DataFrame.from_dict(dosage_dict, orient='index')
    genotype_df.index.name = 'individual'

    genotype_df = genotype_df.loc[genotype_df.index.intersection(samples_list)]
    return genotype_df


def ld_parser(chrom, str_input, egenes_chrom, estrs, gene_lib, samples_list):
    """
    Parses the input VCF file and calculates LD for each eGene in the specified chromosome.
    """
    import pandas as pd
    from cyvcf2 import VCF
    from cpg_utils.config import output_path
    from cpg_utils import to_path

    # Create a VCF reader
    vcf_reader = VCF(str_input['vcf'])

    for gene_name in egenes_chrom['gene_name'].unique():
        if to_path(output_path(f'ld_files/{gene_name}.csv')).exists():
            print(f'LD file for {gene_name} already exists, skipping...')
            continue
        gene_start = max(gene_lib[gene_lib['gene_ids'] == gene_name]['start'].values[0] - 100_000, 0).astype(int)
        gene_end = (gene_lib[gene_lib['gene_ids'] == gene_name]['end'].values[0] + 100_000).astype(int)

        estrs_egene = estrs[estrs['gene_name'] == gene_name]
        df = build_tr_genotype_matrix(vcf_reader, chrom, gene_start, gene_end, estrs_egene, samples_list)
        ld_matrix = df.corr(method='pearson') ** 2
        ld_matrix.to_csv(output_path(f'ld_files/{gene_name}.csv'))



@click.option(
    '--estrs-path',
    type=str,
    required=True,
    help='Path to the eSTRs CSV file.',
    default='gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/cell-type-spec/estrs.csv',
)
@click.option(
    '--tr-vcf-dir',
    type=str,
    required=True,
    help='Directory containing TR VCF files.',
    default='gs://cpg-tenk10k-test/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific',
)
@click.command()
def main(estrs_path, tr_vcf_dir):
    b = get_batch(name='LD matrix runner')
    # load estrs
    estrs = pd.read_csv(estrs_path)
    estrs['variantid'] = estrs['pos'].astype(str) + estrs['motif']

    # load sample list
    bioheart_ids = pd.read_csv(
        'gs://cpg-bioheart-test/tenk10k/str/associatr/final-freeze/input_files/bioheart_n975_sample_covariates.csv'
    )['sample_id']
    tob_ids = pd.read_csv(
        'gs://cpg-bioheart-test/tenk10k/str/associatr/final_freeze/input_files/tob_n950/covariates/6_rna_pcs/CD4_TCM_covariates.csv'
    )['sample_id']
    samples_list = bioheart_ids.to_list() + tob_ids.to_list()
    samples_list = [s.removeprefix('CPG') for s in samples_list]

    # load egene annotations
    gene_lib = pd.read_csv('gs://cpg-tenk10k-test/saige-qtl/300libraries_n1925_adata_raw_var.csv')
    gene_lib = gene_lib[['gene_ids', 'start', 'end']]

    egenes = estrs.drop_duplicates(['chr', 'gene_name'])
    for chrom in [estrs['chr'].unique()]:
        str_vcf_path = f'{tr_vcf_dir}/hail_filtered_{chrom}.vcf.bgz'

        # run LD calculation for each chrom-celltype combination
        ld_job = b.new_python_job(
            f'LD calc for {chrom}',
        )
        str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'})
        egenes_chrom = egenes[egenes['chr'] == chrom]
        ld_job.call(ld_parser, chrom, str_input, egenes_chrom, estrs, gene_lib, samples_list)

    b.run(wait=False)


if __name__ == '__main__':
    main()
