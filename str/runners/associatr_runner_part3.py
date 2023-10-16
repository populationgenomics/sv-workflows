#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member
"""
This script is step 3 of 4 for running associaTR.
It aims to:
- apply locus filters to mergedSTR VCF
- remove "CPG" prefix from CPG IDs in the VCF
- bgzip and tabix the filtered VCF for input into associatr

 analysis-runner --dataset "tob-wgs" \
    --description "Run dumpSTR" \
    --access-level "test" \
    --output-dir "hoptan-str/associatr" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:587e9cf9dc23fe70deb56283d132e37299244209 \
    associatr_runner_part3.py --file-path=gs://cpg-tob-wgs-test/str/expansionhunter/v4-mergeSTR/chr22/mergeSTR_1057_samples_eh.vcf.gz \
    --filter-regions=gs://cpg-tob-wgs-test/hoptan-str/associatr/input_files/segDupRegions/segDupRegions_hg38_sorted.bed.gz


"""
import click

from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch
from cpg_utils import to_path

from cpg_utils.hail_batch import output_path

config = get_config()
TRTOOLS_IMAGE = config['images']['trtools']
BCFTOOLS_IMAGE = config['images']['bcftools']

def file_parser(file_path, ofile_path):
    with to_path(file_path).open('r') as infile, open(ofile_path, 'w') as outfile:
        for line in infile:
            modified_line = line.replace("CPG", "")
            outfile.write(modified_line)

# inputs:
# file-path
@click.option('--file-path', help='gs://... to the output of mergedSTR')
# caller
@click.option(
    '--vcftype',
    help='gangstr or eh',
    type=click.Choice(['eh', 'gangstr'], case_sensitive=True),
    default='eh',
)
# locus level filters
@click.option(
    '--min-locus-call-rate',
    help=' Minimum locus call rate (provided as a proportion eg 0.9)',
    default=0.9,
    type=float,
)
@click.option(
    '--min-locus-het',
    help='Minimum locus heterozygosity (provided as a proportion eg 0.1)',
    default=0.1,
    type=float,
)
@click.option(
    '--min-locus-hwep',
    help='Minimum locus Hardy-Weinberg equilibrium p-value (provided as a proportion eg 0.0001)',
    default=0.0001,
    type=float,
)
@click.option(
    '--filter-regions',
    help='Regions to filter out of the VCF',
    default=None,
    type=str,
)
@click.command()
def main(file_path, vcftype, min_locus_call_rate, min_locus_het, min_locus_hwep, filter_regions):

    b = get_batch()
    merged_str_vcf =b.read_input_group(
        **dict(
            base = file_path,
            tbi = file_path + '.tbi',
        )
    )
    trtools_job = b.new_job(name = f'Run dumpSTR locus level filters')
    trtools_job.storage('20G')
    trtools_job.cpu(4)
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.declare_resource_group(ofile={'vcf': '{root}.vcf',
                                              'loclog.tab': '{root}.loclog.tab',
                                              'samplog.tab':'{root}.samplog.tab'
                                              })

    if filter_regions is not None:
        filter_regions_input = b.read_input_group(
            **dict(
                base = filter_regions,
                tbi = filter_regions + '.tbi',
            )
        )
        trtools_job.command(f' dumpSTR --vcf {merged_str_vcf.base} --out {trtools_job.ofile} --vcftype {vcftype} --min-locus-callrate {min_locus_call_rate} --min-locus-het {min_locus_het} --min-locus-hwep {min_locus_hwep} --filter-regions {filter_regions_input.base} --use-length --drop-filtered')
    else:
        trtools_job.command(f' dumpSTR --vcf {merged_str_vcf.base} --out {trtools_job.ofile} --vcftype {vcftype} --min-locus-callrate {min_locus_call_rate} --min-locus-het {min_locus_het} --min-locus-hwep {min_locus_hwep} --use-length --drop-filtered')

    b.write_output(trtools_job.ofile, output_path(f'input_files/dumpSTR/dumpSTR_filtered'))

    file_parser_job = b.new_python_job(name = f'file_parser')
    file_parser_job.image(config['workflow']['driver_image'])
    file_parser_job.storage('20G')
    file_parser_job.cpu(4)
    file_parser_job.call(file_parser, trtools_job.ofile['vcf'], file_parser_job.ofile)

    bcftools_job = b.new_job(name = f'bgzip and tabix the dumpSTR output VCF')
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.storage('20G')
    bcftools_job.cpu(4)

    bcftools_job.declare_resource_group(
        vcf_output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    bcftools_job.command(
        f"""
    set -ex;
    echo "Compressing";
    bcftools sort {file_parser_job.ofile} | bgzip -c > {bcftools_job.vcf_output['vcf.gz']};

    echo "indexing {bcftools_job.vcf_output['vcf.gz']}";
    tabix -p vcf {bcftools_job.vcf_output['vcf.gz']};
"""
    )
    b.write_output(bcftools_job.vcf_output, output_path(f'input_files/dumpSTR/dumpSTR_filtered'))

    #dumpSTR --vcf mergeSTR_1057_samples_eh.vcf.gz --out filtered_mergeSTR_results --vcftype eh --min-locus-callrate 0.9 --min-locus-het 0.1 --min-locus-hwep 0.0001 --filter-regions segDupRegions_hg38_sorted.bed.gz


    b.run(wait=False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter







