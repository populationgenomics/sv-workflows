#!/usr/bin/env python3
"""
This script concatenates the results of running `coloc_runner.py` (output is per gene) into a single CSV file.

analysis-runner --dataset "bioheart" \
    --description "Parse coloc results" \
    --memory=32G \
    --access-level "test" \
    --output-dir "str/associatr" \
    filter_file.py
"""

import click

import pandas as pd



def main():
    df = pd.read_csv('gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/all_cell_types.tsv', sep = '\t')
    df_subset = df[(df['susie_pip']>=0.8) &(df['finemap_prob']>=0.8)]
    df_subset.to_csv('gs://cpg-bioheart-test-analysis/str/associatr/fine_mapping/susie_finemap/all_cell_types_susie_finemap_0.8.tsv', sep = '\t', index = False)

if __name__ == '__main__':
    main()
