#!/usr/bin/env python3

"""

This script creates a new VCF of just the FM'ed eSTRs that we want to do the GWAS on.

analysis-runner --dataset "bioheart" \
    --description "Create a VCF of just the FM'ed eSTRs" \
    --access-level "test" \
    --memory='4G' \
    --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
    --output-dir "str/associatr/gwas-cell-prop" \
    fm_vcf_maker.py

"""

import pandas as pd
from cyvcf2 import VCF, Writer

import hail as hl

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, init_batch, output_path


def get_fm_estr_vcf(str_input):
    init_batch()

    fm_estr = pd.read_csv('gs://cpg-bioheart-test/str/associatr/gwas-cell-prop/input_files/estrs_lead_filtered.csv')
    fm_estr = fm_estr.drop_duplicates(subset=['chr', 'pos', 'motif_x'])
    fm_estr['locus'] = (
        fm_estr['chr'].astype(str) + ':' + fm_estr['pos'].astype(str) + ':' + fm_estr['motif_x'].astype(str)
    )
    fm_estr_set = set(fm_estr['locus'])
    str_vcf = VCF(str_input['vcf'])

    output_vcf = 'subset.vcf'
    writer = Writer(output_vcf, str_vcf)

    for variant in str_vcf:
        variant_id = variant.CHROM + ':' + str(variant.POS) + ':' + str(variant.INFO.get('RU'))
        if variant_id in fm_estr_set:
            writer.write_record(variant)

    writer.close()
    hl.hadoop_copy(output_vcf, 'gs://cpg-bioheart-test/str/associatr/gwas-cell-prop/input_files/fm_estr.vcf')


def main():
    b = get_batch(name='fm_vcf_maker')

    str_vcf_path = 'gs://cpg-bioheart-test/str/associatr/input_files/vcf/v1/hail_filtered.vcf.bgz'
    str_input = get_batch().read_input_group(**{'vcf': str_vcf_path, 'tbi': str_vcf_path + '.tbi'})

    ld_job = b.new_python_job('get_fm_estr_vcf')
    ld_job.cpu(2)
    ld_job.storage('50G')
    ld_job.call(get_fm_estr_vcf, str_input)
    b.run(wait=False)


if __name__ == '__main__':
    main()
