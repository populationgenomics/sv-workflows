#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member,unnecessary-lambda

import pandas as pd

def main():

    df1 = pd.read_csv('gs://cpg-bioheart-test/str/Trujillo_methylation_eQTLs/hg38_STRs_parsed.tsv', sep = '\t')
    df2 = pd.read_csv('gs://cpg-bioheart-test/str/Trujillo_methylation_eQTLs/hg38_SNVs_parsed.tsv', sep = '\t')

    combined_df = pd.concat([df1, df2], ignore_index=True)

    combined_df.to_csv('gs://cpg-bioheart-test/str/Trujillo_methylation_eQTLs/hg38_STRs_SNVs_parsed.tsv', sep = '\t', index = False)

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter


