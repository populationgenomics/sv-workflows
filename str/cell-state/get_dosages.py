#!/usr/bin/env python3

"""
This script extracts dosages for TRs falling in the cis window of a set of genes.

analysis-runner --dataset tenk10k --access-level test --description "get dosages" \
--output-dir "str/cellstate/input_files" get_dosages.py
"""

import click
import json
from cpg_utils import to_path
import pandas as pd
from cpg_utils.hail_batch import get_batch


def gene_cis_window_file_reader(file_path):
    cis_window = pd.read_csv(file_path, sep='\t', header=None)
    cis_window.columns = ['chrom', 'start', 'end']
    return f'{cis_window["chrom"][0]}:{cis_window["start"][0]}-{cis_window["end"][0]}'


def tr_extract_genotype_matrix(cis_window, vcf):
    """
    Extracts the genotype matrix for tandem repeats from a VCF file.
    (summed repeats)
    """

    import pandas as pd

    region = cis_window
    samples = vcf.samples
    genotype_dict = {sample: {} for sample in samples}

    for variant in vcf(region):
        coord = f"{variant.CHROM}:{variant.POS}:{variant.INFO.get('RU')}"
        repcn_data = variant.format('REPCN')  # array of strings like "10/12"

        if repcn_data is None:
            continue  # skip if field missing

        for i, sample in enumerate(samples):
            val = repcn_data[i]
            try:
                if isinstance(val, bytes):
                    val = val.decode()

                if val in {"./.", ".", ""}:
                    summed = None
                else:
                    allele_1, allele_2 = map(int, val.split('/'))
                    summed = allele_1 + allele_2
            except Exception:
                summed = None

            genotype_dict[sample][coord] = summed

    df = pd.DataFrame.from_dict(genotype_dict, orient='index')
    df.index.name = "sample"
    df.sort_index(axis=1, inplace=True)
    df.reset_index(inplace=True)
    return df


def dosages(chromosome, input_gene_list_dir, cis_window_dir, cell_type, pathway):
    import pandas as pd
    import numpy as np
    from cpg_utils.hail_batch import output_path
    from cyvcf2 import VCF
    from cpg_utils import to_path

    """
    Extracts the dosage files for each gene in the df DataFrame.
    """

    gene_list_path = (
        input_gene_list_dir
        + f'/{pathway}/1_min_pct_cells_expressed/{cell_type}/{chromosome}_{cell_type}_gene_list.json'
    )
    with open(to_path(gene_list_path), 'r') as f:
        gene_list = json.load(f)

    vcf_path = to_path(
        f'gs://cpg-tenk10k-test/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific/hail_filtered_{chromosome}.vcf.bgz'
    )
    vcf_path_index = to_path(
        f'gs://cpg-tenk10k-test/str/associatr/final-freeze/input_files/tr_vcf/v1-chr-specific/hail_filtered_{chromosome}.vcf.bgz.tbi'
    )

    local_vcf_file = 'local.vcf.bgz'
    local_vcf_file_index = 'local.vcf.bgz.tbi'

    vcf_path.copy(local_vcf_file)
    vcf_path_index.copy(local_vcf_file_index)

    vcf_reader = VCF(local_vcf_file)

    for gene in gene_list:
        if to_path(output_path(f'dosages/{chromosome}/{gene}_dosages.csv')).exists():
            print(f"Dosage file for {gene} already exists, skipping.")
            continue

        # need to extract the gene start and end from the cis window file for input into 'region'
        gene_cis_window_file = f'{cis_window_dir}/{pathway}/{cell_type}/{chromosome}/{gene}_100000bp.bed'
        cis_window_region = gene_cis_window_file_reader(gene_cis_window_file)
        df = tr_extract_genotype_matrix(cis_window_region, vcf_reader)
        df.to_csv(output_path(f'dosages/{chromosome}/{gene}_dosages.csv'), index=False)


@click.option('--input-gene-list-dir', default='gs://cpg-tenk10k-test/str/cellstate/input_files/tob/scRNA_gene_lists')
@click.option('--cis-window-dir', default='gs://cpg-tenk10k-test/str/cellstate/input_files/tob/cis_window_files')
@click.option('--cell-type', default='B_naive', help='Cell type to process')
@click.option('--pathway', default='GOBP_MULTI_MULTICELLULAR_ORGANISM_PROCESS_subtype', help='Pathway to process')
@click.option('--chromosome', default='chr22', help='Chromosome to process')
@click.option('--job-storage', help='Storage of the batch job eg 30G', default='8G')
@click.option('--job-memory', help='Memory of the batch job', default='standard')
@click.option('--job-cpu', help='Number of CPUs of Hail batch job', default=2)
@click.command()
def main(input_gene_list_dir, cis_window_dir, cell_type, pathway, chromosome, job_storage, job_memory, job_cpu):
    """
    Run context eQTL analysis for a specific cell type and pathway
    """
    b = get_batch(name='Dosages files prep for Cell state work')
    for chromosome in chromosome.split(','):
        j = b.new_python_job(
            name=f'Extract dosages for {cell_type} {chromosome}',
        )
        j.cpu(job_cpu)
        j.memory(job_memory)
        j.storage(job_storage)
        j.call(
            dosages,
            chromosome,
            input_gene_list_dir,
            cis_window_dir,
            cell_type,
            pathway,
        )
    b.run(wait=False)


if __name__ == '__main__':
    main()
