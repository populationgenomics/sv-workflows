#!/usr/bin/env python3
# pylint: disable=duplicate-code
"""
This script prepares GangSTR/EH VCF files for input into mergeSTR.
Required input: --caller, --input-dir, and external sample IDs
For example:
analysis-runner --access-level test --dataset tob-wgs --description 'tester' --output-dir 'str/5M_run_combined_vcfs/merge_str_prep' merge_str_prep.py --caller=eh --input-dir=gs://cpg-tob-wgs-test/str/5M_run_combined_vcfs CPGtestersorted2

Required packages: sample-metadata, hail, click, os
pip install sample-metadata hail click

"""
import os
import click

from cpg_utils.config import get_config
from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, output_path, reference_path

config = get_config()

REF_FASTA = str(reference_path('broad/ref_fasta'))
BCFTOOLS_IMAGE = config['images']['bcftools']


# input caller
@click.option(
    '--caller',
    help='gangstr or eh',
    type=click.Choice(['eh', 'gangstr'], case_sensitive=True),
)
# input directory
@click.option('--input-dir', help='gs://...')
# input sample ID
@click.argument('internal-wgs-ids', nargs=-1)
# sharded flag
@click.option('--sharded', is_flag=True, help='Assume a sharded catalog was used')
@click.command()
def main(
    caller, input_dir, internal_wgs_ids: list[str], sharded
):  # pylint: disable=missing-function-docstring
    # Initializing Batch
    b = get_batch()

    # Working with CRAM files requires the reference fasta
    ref = b.read_input_group(
        **dict(
            base=REF_FASTA,
            fai=REF_FASTA + '.fai',
            dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
            + '.dict',
        )
    )

    input_vcf_dict = {}

    if sharded:
        input_vcf_dict = {
            id: [str(gspath) for gspath in to_path(f'{input_dir}/{id}').glob('*.vcf')]
            for id in internal_wgs_ids
        }

    else:
        input_vcf_dict = {
            id: [os.path.join(input_dir, f'{id}_{caller}.vcf')]
            for id in internal_wgs_ids
        }

    for id in list(input_vcf_dict.keys()):
        bcftools_job = b.new_job(name=f'{id} {caller} Files prep')
        bcftools_job.image(BCFTOOLS_IMAGE)
        bcftools_job.cpu(4)
        bcftools_job.storage('15G')

        for vcf_file in input_vcf_dict[id]:
            vcf_input = b.read_input(vcf_file)

            if caller == 'eh':
                bcftools_job.declare_resource_group(
                    vcf_sorted={
                        'reheader.vcf.gz': '{root}.reheader.vcf.gz',
                        'reheader.vcf.gz.tbi': '{root}.reheader.vcf.gz.tbi',
                    }
                )
                bcftools_job.command(
                    f"""

                    bcftools reheader -f {ref.fai} {vcf_input} | bcftools sort --temp-dir $BATCH_TMPDIR/ | bgzip -c >  {bcftools_job.vcf_sorted['reheader.vcf.gz']}

                    tabix -f -p vcf {bcftools_job.vcf_sorted['reheader.vcf.gz']}

                    """
                )
                # Output writing
                output_file_name = (vcf_file.split('/')[-1]).split('.')[0]
                output_path_eh = output_path(output_file_name, 'analysis')
                b.write_output(bcftools_job.vcf_sorted, output_path_eh)

            else:
                bcftools_job.declare_resource_group(
                    vcf_sorted={
                        'vcf.gz': '{root}.vcf.gz',
                        'vcf.gz.tbi': '{root}.vcf.gz.tbi',
                    }
                )
                bcftools_job.command(
                    f"""

                    bcftools sort {vcf_input} | bgzip -c  > {bcftools_job.vcf_sorted['vcf.gz']}

                    tabix -p vcf {bcftools_job.vcf_sorted['vcf.gz']}

                    """
                )
                # Output writing
                output_path_gangstr = output_path(f'{id}_gangstr', 'analysis')
                b.write_output(bcftools_job.vcf_sorted, output_path_gangstr)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
