#!/usr/bin/env python3
"""
This script runs FINEMAP.

analysis-runner --dataset "bioheart" --access-level "test" --description "Run finemap" --output-dir "str/associatr/fine_mapping" finemap_runner.py \
    --input-dir=gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/finemap_prep \
    --associatr-dir=gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results \
    --celltypes="ASDC" \
    --chroms="chr22"

"""
import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import get_batch, image_path, output_path


@click.option('--input-file-dir', help='Input file directory')
@click.option('--associatr-dir', help='Directory containing the associatr results')
@click.option('--celltypes')
@click.option('--chroms')
@click.option('--n-causal-snps', help='Number of causal SNPs to be used in the FINEMAP analysis', default=10)
@click.command()
def main(input_file_dir, celltypes, chroms, associatr_dir, n_causal_snps):
    b = get_batch(name=f'Run FINEMAP on {celltypes}:{chroms}')
    for celltype in celltypes.split(','):
        for chrom in chroms.split(','):
            input_files = list(to_path(f'{input_file_dir}/{celltype}/{chrom}').glob('*.z'))
            genes = [str(f).split('/')[-1].split('.')[0] for f in input_files]
            for gene in genes:  # run FINEMAP for each gene
                # load in the associaTR file to extract n_samples tested - used as a FINEMAP command line arg
                associatr_sum_stats = pd.read_csv(
                    f'{associatr_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv', sep='\t',
                )
                n_samples_tested_1 = associatr_sum_stats.loc[0, 'n_samples_tested_1']
                n_samples_tested_2 = associatr_sum_stats.loc[0, 'n_samples_tested_2']
                n_samples_total = n_samples_tested_1 + n_samples_tested_2

                # load in the z and ld files required for FINEMAP
                data_in = b.read_input_group(
                    **{
                        'ld': f'{input_file_dir}/{celltype}/{chrom}/{gene}.ld',
                        'z': f'{input_file_dir}/{celltype}/{chrom}/{gene}.z',
                    },
                )

                finemap_job = b.new_job(name=f'Run finemap on {celltype}:{chrom} {gene}')
                finemap_job.image(image_path('finemap'))
                finemap_job.declare_resource_group(
                    ofile={
                        'data.snp': '{root}.data.snp',
                        'data.config': '{root}.data.config',
                        'data.log': '{root}.data.log',
                        'data.cred': '{root}.data.cred',
                    },
                )
                finemap_job.command(
                    f"""
                            # Create a temp_file that stores the I/O file paths for FINEMAP (as per format required by the tool)
                            temp_file=$(mktemp)

                            # Write the required format to the file
                            echo 'z;ld;snp;config;cred;log;n_samples' > $temp_file
                            echo "{data_in.z};{data_in.ld};{finemap_job.ofile['data.snp']};{finemap_job.ofile['data.config']};data.cred;data.log;{n_samples_total}" >> $temp_file
                            chmod +x $temp_file

                            # Run FINEMAP
                            finemap --sss --in-files $temp_file --log --n-causal-snps {n_causal_snps}

                            # Concatenate all files with data.cred* into a single file
                            cat data.cred* > {finemap_job.ofile['data.cred']}

                            # Cat data.log* into a single file
                            cat data.log* > {finemap_job.ofile['data.log']}

                            """,
                )

                output_path_vcf = output_path(f'finemap/ofiles/{celltype}/{chrom}/{gene}', 'analysis')
                b.write_output(finemap_job.ofile['data.snp'], output_path_vcf + '.snp')
                b.write_output(finemap_job.ofile['data.config'], output_path_vcf + '.config')
                b.write_output(finemap_job.ofile['data.log'], output_path_vcf + '.log')
                b.write_output(finemap_job.ofile['data.cred'], output_path_vcf + '.cred')

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
