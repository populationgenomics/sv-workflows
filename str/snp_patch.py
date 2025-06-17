#!/usr/bin/env python3

"""

analysis-runner --dataset "bioheart" --access-level "test" --description "Patch SNP files" --output-dir "tenk10k/str/associatr/final_freeze/common_variant_snps/tob_n950/results/v2-patch" snp_patch.py \
--chromosomes=chr21 --cell-types=CD4_CTL
"""

import click


from cpg_utils.hail_batch import get_batch, output_path
from cpg_utils import to_path


def correct_csv_format(file_path, cell, chrom):
    from cpg_utils import to_path
    from cpg_utils.hail_batch import output_path
    import pandas as pd

    try:
        df = pd.read_csv(file_path, sep='\t')

    except pd.errors.ParserError:
        corrected_lines = []
        with open(to_path(file_path), "r") as f:
            for i, line in enumerate(f):
                fields = line.strip().split("\t")
                if i == 0:
                    # Assume first line is header
                    header = fields
                    corrected_lines.append(fields)
                elif len(fields) == 13:
                    corrected_lines.append(fields)
                elif len(fields) > 13 and len(fields) >= 18:
                    # Malformed line: ignore first 5 fields, keep the rest to make up 13 total
                    fixed_fields = fields[5:]
                    corrected_lines.append(fixed_fields)

        # Convert to DataFrame using the fixed data
        df = pd.DataFrame(corrected_lines[1:], columns=corrected_lines[0])
        df.to_csv(
            output_path(f'{cell}/{chrom}/{to_path(file_path).name}', 'analysis'),
            sep='\t',
            index=False,
        )


@click.option(
    "--input-dir",
    help="Directory containing the input files",
    default='gs://cpg-bioheart-test-analysis/tenk10k/str/associatr/final_freeze/common_variant_snps/tob_n950/results/v1',
)
@click.option('--cell-types', help='Cell types to process', default='CD4_TCM,CD4_Naive,CD4_TEM,CD4_CTL,CD4_Proliferating,NK,NK_CD56bright,NK_Proliferating,CD8_TEM,CD8_TCM,CD8_Proliferating,CD8_Naive,Treg,B_naive,B_memory,B_intermediate,Plasmablast,CD14_Mono,CD16_Mono,cDC1,cDC2,pDC,dnT,gdT,MAIT,ASDC,HSPC,ILC')
@click.option('--chromosomes', help='Chromosomes to process', default='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22')
@click.command()
def main(input_dir, cell_types, chromosomes):
    b = get_batch(name='SNP file patch')
    for cell_type in cell_types.split(','):
        for chromosome in chromosomes.split(','):
            files = list(to_path(f'{input_dir}/{cell_type}/{chromosome}').glob('*.tsv'))
            for file in files:
                if to_path(output_path(f'{cell_type}/{chromosome}/{to_path(file.name)}', 'analysis')).exists():
                    print(f"File {file} already processed, skipping.")
                    continue
                job = b.new_python_job(f'Patch SNPs for {file.name}')
                job.cpu(0.25)
                job.call(correct_csv_format, str(file), cell_type, chromosome)
    b.run(wait=False)


if __name__ == "__main__":
    main()
