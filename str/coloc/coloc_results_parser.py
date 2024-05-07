#!/usr/bin/env python3
"""
analysis-runner --dataset "bioheart" \
    --memory='16G'
    --description "Parse coloc results" \
    --access-level "test" \
    --output-dir "str/associatr/coloc/CD4_TCM" \
    coloc_results_parser.py


"""
from cpg_utils import to_path
import pandas as pd

def main():
    files = list(to_path(f'gs://cpg-bioheart-test/str/associatr/coloc/CD4_TCM').glob('*.tsv'))
    # List to store DataFrames
    dfs = []

    # Iterate over each file path
    for file_path in files:
        # Read file into a DataFrame
        df = pd.read_csv(file_path, sep = '\t')
        # Append DataFrame to the list
        dfs.append(df)

    # Concatenate DataFrames row-wise
    result_df = pd.concat(dfs, ignore_index=True)
    # Write the result to a CSV file
    result_df.to_csv('gs://cpg-bioheart-test/str/associatr/coloc/CD4_TCM/gene_summary_result.csv', index=False)

if __name__ == '__main__':
    main()
