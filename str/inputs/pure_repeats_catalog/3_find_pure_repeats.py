"""
This script identifies pure repeat loci (motif 2-6bp) in the Illumina catalog, with no
tolerance for indels/SNPs. It then creates left and right flanking sequences of pure
repeat loci, which will be used as input in `pure_repeats_catalog.py` to see if there
are additional repeats in the flanks.

Required packages: str_analysis
python3 -m pip install --upgrade str_analysis
# TODO this is best contained in a requirements.txt file at the project root
"""


from str_analysis.utils.find_repeat_unit import find_repeat_unit_kmer


INPUT_FILE = 'intermediate_files/Illumina_catalog_sequences.fasta.txt'
EXCLUDED_LOCUS_FILE = 'intermediate_files/excluded_loci_catalog.txt'
PURE_CATALOG_OUTPUT = "intermediate_files/pure_repeat_catalog_not_final.txt"
CATALOG_ONE_BED = "intermediate_files/catalog_with_flanks_one_motif_length.bed"
CATALOG_TWO_BED = "intermediate_files/catalog_with_flanks_two_motif_length.bed"


catalog_dict = {}
excluded_dict = {}

# manual removal of 2 loci that overlap each other
# (representing a compound repeat structure)
manually_excluded_loci = {'chr8:49170487-49170547', 'chr8:49170544-49170552'}

# region: read input file
with open(INPUT_FILE) as handle:
    for line in handle:

        # find an identifier line - use it to find the paired seq
        if line.startswith('>'):
            # remove the > and trailing whitespace
            key = line.removeprefix('>').rstrip()
            # pull out the next line, read both together
            value = next(handle).rstrip()

            # get the kmer value for this sequence
            repeat_unit, repeat_count = find_repeat_unit_kmer(value)

            # homopolymer
            # impure repeat sequence
            # motifs greater than 6bp (STR definition is 2-6bp motif)
            # manually excluded sites
            if ((repeat_count == 1) or (len(repeat_unit) == 1) or (len(repeat_unit) > 6)) or key in manually_excluded_loci:
                excluded_dict[key] = (repeat_unit, repeat_count)
            else:
                catalog_dict[key] = (repeat_unit, repeat_count)
# endregion

# region: write out the excluded loci
with open(EXCLUDED_LOCUS_FILE, 'w') as handle:
    for key, (motif, repeat_count) in excluded_dict.items():
        handle.write(f'{key} {motif} {repeat_count}\n')
# endregion

print(f'Number of pure repeats: {len(catalog_dict.keys())}')  # expected 164847 pure repeat loci

# region: write out the pure repeat loci to an output file
with open(PURE_CATALOG_OUTPUT, 'w') as handle:
    for key, (motif, repeat_count) in catalog_dict.items():
        handle.write(f'{key} {motif} {repeat_count}\n')
# endregion


# region: write pure loci with flanking sequences
# opens two files, iterates over pure loci once, writes out both variations
with open(CATALOG_ONE_BED, 'w') as one_handle:
    with open(CATALOG_TWO_BED, 'w') as two_handle:
        for key, (motif, repeat_count) in catalog_dict.items():
            chrom, numerical = key.split(':')[0]

            # split the start and end coordinates, apply int() to both
            (start, end) = map(int, numerical.split('-'))
            motif_len = len(motif)

            # write out the original sequence region
            one_handle.write(f'{chrom}\t{start}\t{end}\n')
            two_handle.write(f'{chrom}\t{start}\t{end}\n')

            # write out the sequence region with the right flank
            one_handle.write(f'{chrom}\t{start}\t{end+motif_len}\n')
            two_handle.write(f'{chrom}\t{start}\t{end+(2*motif_len)}\n')

            # write out the sequence region with the left flank
            one_handle.write(f'{chrom}\t{start-motif_len}\t{end}\n')
            two_handle.write(f'{chrom}\t{start-(motif_len*2)}\t{end}\n')
# endregion
