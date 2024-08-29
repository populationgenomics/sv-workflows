#!/usr/bin/env python3


""""
Raw methylation calls from PacBio data (using pb-CpG-tools) are provided in a BED per-sample format.

This script concatenates each per-sample BED file into a single BED file (chromosome-specific chr1-22).
Only CpG sites called in all samples are retained (inner join merge).
Only the 'mod_score' parameter is extracted; all other methylation-related columns are discarded.

analysis-runner --dataset "bioheart" --access-level "test" --description "Concatenate methylation BED files" --output-dir "str/pacbio-methylation/combined_bed" methylation_bed_parser.py


"""

import click


from cpg_utils.hail_batch import get_batch



def concatenator(input_methylation_dir, chrom_num):
    import pandas as pd
    from cpg_utils import to_path

    from cpg_utils.hail_batch import output_path
    methylation_files = list(to_path(f'{input_methylation_dir}').glob('*.combined.bed'))
    chrom = f'chr{chrom_num}'
    print(f'Processing {chrom}...')
    master_df = pd.DataFrame()
    for file in methylation_files:
        file_name = str(file)
        sample = file_name.split('/')[-1].split('.')[0]
        df = pd.read_csv(
            file, sep='\t', usecols=[0, 1, 3], names=['chrom', 'start', f'{sample}'],
        )  # col3 corresponds to mod_score
        df = df[df['chrom'] == chrom]
        if master_df.empty:
            master_df = df
        else:
            master_df = master_df.merge(df, on=['chrom', 'start'], how='inner')
    output_gcs = output_path(f'methylation_combined_{chrom}.bed')
    master_df.to_csv(output_gcs, sep='\t', index=False, header=True)

@click.option(
    '--input-methylation-dir',
    help='GCS path to the directory containing the methylation BED files',
    default='gs://cpg-bioheart-test/str/pacbio-methylation',

)
@click.command()
def main(input_methylation_dir):
    b = get_batch(name='Methylation bed parser')
    for chrom_num in range(1, 23):
        methylation_parser_job = b.new_python_job(f'Methylation parser for chr{chrom_num}')
        methylation_parser_job.cpu(2)
        methylation_parser_job.memory('20G')
        methylation_parser_job.call(concatenator, input_methylation_dir, chrom_num)
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
