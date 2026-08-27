#!/usr/bin/env python3
# pylint: disable=no-value-for-parameter
"""
Extract the raw p-values from the output of associatr_runner_sv.py into one TSV per cell type,
for downstream QQ plots / lambda_GC. SV counterpart of str/associatr/helper/raw_pval_extractor.py.

Differences from the STR version, all forced by the SV output format:

1. Columns are located by NAME, not by position. The STR script hardcodes indices (pval at
   iloc[:, 5]); associaTR names the p-value column after the phenotype
   (`p_{celltype}_{chrom}_{gene}`), so it is found by the `p_` prefix instead. Same convention as
   run_gene_level_pvals_meta_sv.py.

2. Each row carries SV-appropriate locus attributes. An SV is an interval, so POS alone does not
   identify it: associaTR stores `<REF>-<ALT>-END:<end>-<VARID>` in the `motif` column, and that is
   unpacked into `end`, `svtype` and `varid` columns. `varid` is what joins back to the VCF and to
   the federated pairing table used by meta_runner_sv.py.

3. Loci associaTR did not test are dropped (`locus_filtered != False`, or a NaN p-value). The STR
   script wrote those rows out with a literal 'NA' p-value, which a QQ plot cannot use. Counts of
   what was dropped are printed per cell type.

4. Gene names come from stripping the `_{cis_window_size}bp` suffix, not from splitting the
   filename on the first '_'.

Input layout is what associatr_runner_sv.py writes:
    {input-dir}/{cell_type}/{chrom}/{gene}_{cis_window_size}bp.tsv

analysis-runner --dataset "tenk10k" --description "SV raw pval extractor" --access-level "test" \
    --output-dir "tenk10k-sv/sv/sensitivity_analysis/tob_n935/7pc/v1" \
    python3 raw_pval_extractor_sv.py \
    --input-dir=gs://cpg-tenk10k-test-analysis/tenk10k-sv/sv/sensitivity_analysis/tob_n935/7pc/results/v1 \
    --cell-types=CD4_TCM,CD4_TCM_permuted \
    --chromosomes=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
"""

import re
from concurrent.futures import ThreadPoolExecutor

import click
import pandas as pd

from cpg_utils import to_path
from cpg_utils.hail_batch import output_path

# associaTR RU/motif field is '<REF>-<ALT>-END:<end>-<VARID>', e.g.
# 'N-<INS>-END:168956831-ALL_BATCHES.CHR1.FINAL_CLEANUP_INS_CHR1_2102'
_MOTIF_RE = re.compile(r'^[^-]*-<?(?P<svtype>[^>-]+)>?-END:(?P<end>\d+)-(?P<varid>.+)$')

# trailing '_1000000bp' on the per-gene result files
_GENE_RE = re.compile(r'_\d+bp$')

OUTPUT_COLUMNS = ('chrom', 'pos', 'end', 'svtype', 'varid', 'gene', 'pval')


def _gene_name(gene_file: str) -> str:
    """'.../ENSG00000000457_1000000bp.tsv' -> 'ENSG00000000457'."""
    return _GENE_RE.sub('', str(gene_file).split('/')[-1].removesuffix('.tsv'))


def _unpack_motif(motif):
    """'N-<DUP>-END:168957979-ALL_BATCHES...DUP_CHR1_912' -> ('DUP', '168957979', 'ALL_BATCHES...')."""
    match = _MOTIF_RE.match(str(motif))
    if not match:
        return 'NA', 'NA', 'NA'
    return match.group('svtype'), match.group('end'), match.group('varid').strip()


def read_gene_file(gene_file: str):
    """Read one associaTR gene TSV; return (rows, n_dropped).

    rows is a list of OUTPUT_COLUMNS-ordered tuples for the loci that were actually tested.
    """
    try:
        with to_path(gene_file).open() as handle:
            gene_results = pd.read_csv(handle, sep='\t', dtype={'locus_filtered': str})
    except (FileNotFoundError, OSError) as e:
        print(f'!! could not read {gene_file}: {e}')
        return [], 0
    if gene_results.empty:
        return [], 0

    # the p-value column is named after the phenotype (p_{celltype}_{chrom}_{gene})
    pval_columns = [c for c in gene_results.columns if c.startswith('p_')]
    if len(pval_columns) != 1:
        raise ValueError(f'expected exactly one p_* column in {gene_file}, found {pval_columns}')

    n_total = len(gene_results)
    # associaTR flags loci it could not test, and reports their p-value as NA
    gene_results = gene_results[gene_results['locus_filtered'] == 'False']
    gene_results = gene_results[gene_results[pval_columns[0]].notna()]

    gene = _gene_name(gene_file)
    rows = []
    for chrom, pos, motif, pval in zip(
        gene_results['chrom'],
        gene_results['pos'],
        gene_results['motif'],
        gene_results[pval_columns[0]],
    ):
        svtype, end, varid = _unpack_motif(motif)
        rows.append((str(chrom), str(pos), end, svtype, varid, gene, str(pval)))
    return rows, n_total - len(rows)


@click.option('--input-dir', required=True, help='GCS path to the results dir from associatr_runner_sv.py')
@click.option('--cell-types', required=True, help='Name of the cell type, comma separated if multiple')
@click.option('--chromosomes', required=True, help='Chromosome eg chr1 or 1, comma separated if multiple')
@click.option('--threads', type=int, default=16, help='concurrent GCS reads')
@click.command()
def main(input_dir, cell_types, chromosomes, threads):
    """Extract the raw p-values from the SV associaTR results into one TSV per cell type."""
    for cell_type in cell_types.split(','):
        gcs_output = output_path(f'raw_pval_extractor/{cell_type}_gene_tests_raw_pvals.tsv', 'analysis')
        n_written = 0
        n_dropped = 0
        n_genes = 0

        with to_path(gcs_output).open('w') as f:
            f.write('\t'.join(OUTPUT_COLUMNS) + '\n')
            for chromosome in chromosomes.split(','):
                # accept either 'chr1' or '1'; associatr_runner_sv.py writes 'chr1'-style directories
                chrom = chromosome if chromosome.startswith('chr') else f'chr{chromosome}'
                gene_files = sorted(map(str, to_path(f'{input_dir}/{cell_type}/{chrom}').glob('*.tsv')))
                print(f'{cell_type} {chrom}: {len(gene_files):,} gene files')
                n_genes += len(gene_files)

                # reads dominate the runtime; executor.map keeps the output ordered by gene
                with ThreadPoolExecutor(max_workers=threads) as executor:
                    for rows, dropped in executor.map(read_gene_file, gene_files):
                        for row in rows:
                            f.write('\t'.join(row) + '\n')
                        n_written += len(rows)
                        n_dropped += dropped

        print(
            f'{cell_type}: wrote {n_written:,} tested gene-SV pairs from {n_genes:,} genes '
            f'({n_dropped:,} untested loci dropped) -> {gcs_output}',
        )


if __name__ == '__main__':
    main()
