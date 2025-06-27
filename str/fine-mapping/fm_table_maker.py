#!/usr/bin/env python3


"""

This scripts pulls out the coeff_meta_fixed and pval_meta_fixed for each fm TR

analysis-runner --dataset tenk10k --access-level test --output-dir potato --description "make fm table" python3 fm_table_maker.py

"""

import pandas as pd

estrs_fm = pd.read_csv('gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/finemapped_etrs_v2.csv')
estrs_fm[['chr', 'pos', 'motif']] = estrs_fm['variant_id'].str.split(':', expand=True)

# Create empty list to collect each eqtl_row
rows_list = []

for row in estrs_fm.itertuples():
    cell_type = row.cell_type
    gene = row.gene_name
    chrom = row.chr
    pos = int(row.pos)
    motif = row.motif

    eqtl_file_path =f'gs://cpg-tenk10k-test-analysis/str/associatr/final_freeze/tob_n950_and_bioheart_n975/trs_snps/rm_str_indels_dup_strs_v3/{cell_type}/{chrom}/{gene}_100000bp_meta_results.tsv'
    eqtl = pd.read_csv(eqtl_file_path, sep='\t')

    eqtl = eqtl[(eqtl['motif'] == motif) & (eqtl['pos'] == pos)]
    if eqtl.empty:
        print(f'No eQTL data for {gene} in {cell_type} at {chrom}:{pos}_{motif}')
        continue

    eqtl['pval_meta_fixed'] = eqtl['pval_meta_fixed'].astype(float)
    eqtl_row = eqtl.sort_values('pval_meta_fixed').iloc[0]
    eqtl_row['gene_name']= gene
    eqtl_row['cell_type'] = cell_type

    # Append as DataFrame row (convert Series to DataFrame first)
    rows_list.append(eqtl_row.to_frame().T)

# Concatenate all rows
final_df = pd.concat(rows_list, ignore_index=True)

# Write to CSV
final_df.to_csv('gs://cpg-tenk10k-test/str/associatr/final_freeze/meta_fixed/finemapped_etrs_v2_rows.csv', index=False)


