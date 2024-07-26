#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,too-many-arguments
"""
This script runs associaTR (GWAS style)

analysis-runner --dataset "bioheart" --description "run associatr" \
    --access-level "test" \
    --output-dir "str/associatr/gwas-cell-prop/tob_n1055" \
    associatr_runner.py \
    --pheno-cov-dir "gs://cpg-bioheart-test/str/associatr/gwas-cell-prop/tob_n1055/input_files/pheno_cov_numpy"
"""

import click
import numpy as np
import pandas as pd

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, output_path


@click.option('--vcf-file-path', default='gs://cpg-bioheart-test/str/associatr/gwas-cell-prop/input_files/fm_estr.vcf')
@click.option('--pheno-cov-dir', help='Directory containing the phenotype and covariate files')
@click.command()
def main(vcf_file_path, pheno_cov_dir):
    b = get_batch(name='Run associatr')

    for macrogroup in [
        'B',
        'CD4T',
        'CD8T',
        'Cycling',
        'DC',
        'MAIT',
        'Mono',
        'NK',
        'Platelet',
        'Treg',
        'gdT',
        'nan',
        'Eryth',
    ]:
        variant_vcf = b.read_input_group(**{'vcf': vcf_file_path})

        # run associatr for the macrogroup proportion
        associatr_job = b.new_job(name=f'associatr_{macrogroup}')
        gene_pheno_cov = b.read_input(f'{pheno_cov_dir}/{macrogroup}_pheno_cov.npy')
        associatr_job.image(get_config()['images']['trtools'])
        associatr_job.declare_resource_group(association_results={'tsv': '{root}.tsv'})
        associatr_job.command(
            f" associaTR {associatr_job.association_results['tsv']} {variant_vcf.vcf} {macrogroup}_proportion {gene_pheno_cov} --vcftype=eh",
        )
        b.write_output(
            associatr_job.association_results,
            output_path(
                f'results/{macrogroup}_fmed_strs',
                'analysis',
            ),
        )
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
