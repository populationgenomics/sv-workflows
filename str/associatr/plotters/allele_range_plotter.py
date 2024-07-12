#!/usr/bin/env python3
"""
This script plots a QQ plot of observed vs expected -log10(p-values) for each cell type.

analysis-runner --dataset "bioheart" --description "plot qq plot" --access-level "test" \
    --output-dir "str/associatr/tob_n1055_and_bioheart_n990" --memory=16G \
    allele_range_plotter.py


"""
import matplotlib.pyplot as plt
import pandas as pd

import hail as hl

from cpg_utils.hail_batch import init_batch


def main():
    init_batch()
    df = pd.read_csv(
        'gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/annotated_mt/v2/alleles_minus_mode_ht.tsv.bgz',
        compression='gzip',
        sep='\t',
    )
    proportion = ((df['alleles_minus_mode'] > 10) | (df['alleles_minus_mode'] < -10)).mean()
    print(f"Proportion of rows where 'alleles_minus_mode' is >10 or <-10: {proportion}")

    proportion_2 = ((df['alleles_minus_mode'] > 20) | (df['alleles_minus_mode'] < -20)).mean()
    print(f"Proportion of rows where 'alleles_minus_mode' is >10 or <-10: {proportion_2}")
    # Create QQ plot
    plt.figure(figsize=(12, 6))
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.hist(
        df['alleles_minus_mode'],
        bins=30,
        alpha=0.5,
        color='darkblue',
        edgecolor='black',
        density=True,
        range=(-20, 20),
    )

    # Add titles and labels
    plt.xlabel('Allele relative to mode allele', fontsize=14)
    plt.ylabel('Proportion of calls', fontsize=14)
    plt.xticks(fontsize=12)  # Rotate x labels if needed for readability
    plt.yticks(fontsize=12)

    gcs_output_path = 'gs://cpg-bioheart-test/str/wgs_genotyping/polymorphic_run_n2045/plots/alleles_range.png'
    fig.tight_layout()
    fig.savefig('qqplot.png')
    hl.hadoop_copy('qqplot.png', gcs_output_path)


if __name__ == '__main__':
    main()
