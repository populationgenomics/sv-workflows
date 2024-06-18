#!/usr/bin/env python3
"""
This script runs FINEMAP.

analysis-runner --dataset "bioheart" --access-level "test" --description "Run finemap" --output-dir "str/associatr/fine_mapping" \
     --image "australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:d4922e3062565ff160ac2ed62dcdf2fba576b75a-hail-8f6797b033d2e102575c40166cf0c977e91f834e" \
     finemap_runner.py \
    --input-file-dir=gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/finemap_prep \
    --associatr-dir=gs://cpg-bioheart-test-analysis/str/associatr/snps_and_strs/rm_str_indels_dup_strs/v2-whole-copies-only/tob_n1055_and_bioheart_n990/meta_results \
    --celltypes="ASDC" \
    --chroms="chr1,chr2,chr3,chr4,chr5,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22" \

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
@click.option('--job-cpu', help='Number of CPUs to use for each job', default=0.25)
@click.command()
def main(input_file_dir, celltypes, chroms, associatr_dir, n_causal_snps, job_cpu):
    b = get_batch(name=f'Run FINEMAP on {celltypes}:{chroms}')
    for celltype in celltypes.split(','):
        for chrom in chroms.split(','):
            input_files = list(to_path(f'{input_file_dir}/{celltype}/{chrom}').glob('*.z'))
            genes = [str(f).split('/')[-1].split('.')[0] for f in input_files]
            for gene in genes:  # run FINEMAP for each gene
                if to_path(output_path(f'finemap/ofiles/{celltype}/{chrom}/{gene}.snp', 'analysis')).exists():
                    continue
                n_casual_snps_gene = n_causal_snps
                # load in the associaTR file to extract n_samples tested - used as a FINEMAP command line arg
                associatr_sum_stats = pd.read_csv(
                    f'{associatr_dir}/{celltype}/{chrom}/{gene}_100000bp_meta_results.tsv',
                    sep='\t',
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

                data_z = pd.read_csv(f'{input_file_dir}/{celltype}/{chrom}/{gene}.z', sep=' ')
                num_rows = data_z.shape[0]
                if num_rows < n_causal_snps:
                    n_casual_snps_gene = num_rows
                if num_rows == 1:
                    # read in z file
                    data_z = pd.read_csv(f'{input_file_dir}/{celltype}/{chrom}/{gene}.z', sep=' ')

                    # Create a new DataFrame with the desired headers
                    output_data = pd.DataFrame(
                        columns=[
                            'index',
                            'rsid',
                            'chromosome',
                            'position',
                            'allele1',
                            'allele2',
                            'maf',
                            'beta',
                            'se',
                            'z',
                            'prob',
                            'log10bf',
                            'mean',
                            'sd',
                            'mean_incl',
                            'sd_incl',
                        ],
                    )

                    # Create a new row with NA values
                    new_row = pd.Series(
                        {
                            'index': 'NA',
                            'rsid': data_z['rsid'][0],
                            'chromosome': data_z['chromosome'][0],
                            'position': data_z['position'][0],
                            'allele1': 'NA',
                            'allele2': 'NA',
                            'maf': 'NA',
                            'beta': data_z['beta'][0],
                            'se': data_z['se'][0],
                            'z': 'NA',
                            'prob': 1,
                            'log10bf': 'NA',
                            'mean': 'NA',
                            'sd': 'NA',
                            'mean_incl': 'NA',
                            'sd_incl': 'NA',
                        },
                    )

                    # Append the new row to the output DataFrame
                    output_data = output_data.append(new_row, ignore_index=True)
                    output_data.to_csv(
                        output_path(f'finemap/ofiles/{celltype}/{chrom}/{gene}.snp', 'analysis'), sep=' ', index=False,
                    )
                    continue

                finemap_job = b.new_job(name=f'Run finemap on {celltype}:{chrom} {gene}')
                finemap_job.image(image_path('finemap'))
                finemap_job.cpu(job_cpu)
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
                            finemap --sss --in-files $temp_file --log --n-causal-snps {n_casual_snps_gene}

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
