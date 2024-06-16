#!/usr/bin/env python3
"""
This script runs FINEMAP.

analysis-runner --dataset "bioheart" --access-level "test" --description "Run finemap" --output-dir "str/finemapping" finemap_runner.py

"""
from cpg_utils.hail_batch import get_batch, output_path
import click
from cpg_utils import to_path


@click.option('--input-file-dir', help='Input file directory')
@click.option('--celltypes')
@click.option('--chroms')
@click.command()

def main(input_file_dir, celltypes, chroms):
    b = get_batch(name = f'Run FINEMAP on {celltypes}:{chroms}')
    for celltype in celltypes.split(','):
        for chrom in chroms.split(','):
            input_files =list(to_path(f'{input_file_dir}/{celltype}/{chrom}').glob('*.z'))
            genes = [str(f).split('/')[-1].split('.')[0] for f in input_files]
            for gene in genes: #run FINEMAP for each gene
                data_in = b.read_input_group(**{
                        'ld':f'{input_file_dir}/{celltype}/{chrom}/{gene}.ld',
                        'z' : f'{input_file_dir}/{celltype}/{chrom}/{gene}.z'})



                finemap_job = b.new_job(name=f'Run finemap on {celltype}:{chrom} {gene}')
                finemap_job.image('australia-southeast1-docker.pkg.dev/cpg-common/images-dev/finemap:1.4.2')
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
                            # Create a temporary file
                            temp_file=$(mktemp)

                            # Write the required format to the file
                            echo 'z;ld;snp;config;cred;log;n_samples' > $temp_file
                            echo "{data_in.z};{data_in.ld};{finemap_job.ofile['data.snp']};{finemap_job.ofile['data.config']};data.cred;data.log;5363" >> $temp_file
                            chmod +x $temp_file
                            finemap --sss --in-files $temp_file --dataset 1 --log

                            # Concatenate all files with data.cred* into a single file
                            cat data.cred* > {finemap_job.ofile['data.cred']}

                            # Cat data.log* into a single file
                            cat data.log* > {finemap_job.ofile['data.log']}




                            """,
                            )


                output_path_vcf = output_path(f'finemap/example/ofiles')
                b.write_output(finemap_job.ofile['data.snp'], output_path_vcf+'data.snp')
                b.write_output(finemap_job.ofile['data.config'], output_path_vcf+'data.config')
                b.write_output(finemap_job.ofile['data.log'], output_path_vcf+'data.log')
                b.write_output(finemap_job.ofile['data.cred'], output_path_vcf+'data.cred')
    b.run(wait=False)
if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter