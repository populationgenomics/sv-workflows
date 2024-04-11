"""
This script identifies if there are additional repeats in the flanks of pure repeat
loci, and modifies the coordinates to include them as necessary.

Outputs a bed file containing the pure repeat loci coordinates.

Required packages: str_analysis
python3 -m pip install --upgrade str_analysis
"""

import json

from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence
from utils import break_coordinate_string, read_lines_from_file

INPUT_FASTA = 'intermediate_files/catalog_with_flanks_one_motif_length.fasta'
PURE_NOT_FINAL_JSON = 'intermediate_files/pure_repeat_catalog_not_final.json'
LEFT_FLANK_ONE = 'intermediate_files/left_flank_one_motif_added_repeat_loci.txt'
TWO_MOTIF_FLANKS = 'intermediate_files/catalog_with_flanks_two_motif_lengths.fasta'
PURE_BED = 'intermediate_files/pure_repeats_loci.bed'


def read_fasta_content_from_6_line_groups(filename: str) -> dict[str, tuple[str, str]]:
    """
    take a file name, open for reading
    - read contents in groups of 6 lines
    - first is key
    - 4th & 6th are returned as values in a tuple
    :param filename:
    :return:
    """
    # create something to return the data
    return_dict = {}

    with open(filename, encoding='utf-8') as read_handle:
        while True:
            line_group = read_lines_from_file(read_handle, 6)
            if not line_group:
                break

            # groom the strings
            line_group = [line.rstrip() for line in line_group]

            # find an identifier line
            if not line_group[0].startswith('>'):
                raise ValueError('Expected identifier line')

            # remove the > and trailing whitespace
            key = line_group[0].removeprefix('>')

            # get the sequence with right and left flanking sequence
            value = (line_group[3], line_group[5])

            # populate values into the dictionary
            return_dict[key] = value

    return return_dict


# region: read input fasta file
flank_catalog_one_motif_length_dict = read_fasta_content_from_6_line_groups(INPUT_FASTA)
# endregion

# region: load pure catalog dictionary
with open(PURE_NOT_FINAL_JSON, encoding='utf-8') as handle:
    pure_catalog_dict = json.load(handle)
# endregion

# region: look for additional repeats
left_flank_one_motif_added_repeat_loci = []
right_flank_one_motif_added_repeat_loci = []
both_flanks_one_motif_added_repeat_loci = []

for locus, (right_flank, left_flank) in flank_catalog_one_motif_length_dict.items():
    motif, num_repeats = pure_catalog_dict[locus]

    # returns the num_repeats of the longest stretch of repeats with the
    # motif. Sequence input must start with the motif, else returns 0.
    right_repeats = extend_repeat_into_sequence(motif, right_flank)[0]
    left_repeats = extend_repeat_into_sequence(motif, left_flank)[0]

    if left_repeats > int(num_repeats) and right_repeats > int(num_repeats):
        both_flanks_one_motif_added_repeat_loci.append(locus)
    elif left_repeats > int(num_repeats):
        left_flank_one_motif_added_repeat_loci.append(locus)
    elif right_repeats > int(num_repeats):
        right_flank_one_motif_added_repeat_loci.append(locus)

print(f'left flank one added: {len(left_flank_one_motif_added_repeat_loci)}')  # 31
print(f'right flank one added: {len(right_flank_one_motif_added_repeat_loci)}')  # 0
print(f'both flanks one added: {len(both_flanks_one_motif_added_repeat_loci)}')  # 0
# endregion

# region: write left flank added loci
# right_flank_added_repeat_loci and both_flanks_added_repeat_loci are empty
with open(LEFT_FLANK_ONE, 'w', encoding='utf-8') as handle:
    for locus in left_flank_one_motif_added_repeat_loci:
        handle.write(f'{locus}\n')
# endregion

# Because there were some additional repeats detected in the left flanks, will
# consider L/R flanks that are 2 motifs long now.

# File loading of catalog with L/R flanks that are two motif lengths long.
flank_catalog_two_motif_lengths_dict = read_fasta_content_from_6_line_groups(TWO_MOTIF_FLANKS)
print(f'flank_catalog_two_motif_lengths_dict (164847): {len(flank_catalog_two_motif_lengths_dict)}')

# Look into L flank for additional repeats
left_flank_two_motif_added_repeat_loci = []

for locus, (_right_flank, left_flank) in flank_catalog_two_motif_lengths_dict.items():
    motif, num_repeats = pure_catalog_dict[locus]
    left_repeats = extend_repeat_into_sequence(motif, left_flank)[0]

    # just checking L flank now because the R flank yielded
    # no additional repeats in prev. iteration
    if left_repeats > int(num_repeats):
        left_flank_two_motif_added_repeat_loci.append(locus)

# 0, suggesting that all the repeats in the left flanks are only one copy
print(len(left_flank_two_motif_added_repeat_loci))

# region: Edit pure repeats catalog to append repeats in L flank
revised_pure_catalog_dict = {}

for locus, (motif, repeat_count) in pure_catalog_dict.items():
    if locus in left_flank_one_motif_added_repeat_loci:
        motif_length = len(motif)
        chrom, start, end = break_coordinate_string(locus)
        revised_locus = f'{chrom}:{start-motif_length}-{end}'
        revised_repeat_count = int(repeat_count) + 1
        revised_pure_catalog_dict[revised_locus] = (motif, revised_repeat_count)
    else:
        revised_pure_catalog_dict[locus] = (motif, repeat_count)
print(len(revised_pure_catalog_dict))  # 164847
# endregion

# region: testing out revised pure catalog:
# (no harm keeping these in as asserts, no manual review needed)
assert pure_catalog_dict['chr10:10022442-10022450'] == ['CTGC', 2]
assert revised_pure_catalog_dict['chr10:10022438-10022450'] == ('CTGC', 3)
assert pure_catalog_dict['chr10:129364172-129364184'] == ['AATA', 3]
assert revised_pure_catalog_dict['chr10:129364168-129364184'] == ('AATA', 4)
# endregion

# Write out revised pure repeat catalog as a BED file in the format
# chr \t start_coordinate \t end_coordinate \t motif_length \t motif \t num_repeats
with open(PURE_BED, 'w', encoding='utf-8') as handle:
    for p, (motif, repeat_count) in revised_pure_catalog_dict.items():
        motif_length = len(motif)
        chrom, start, end = break_coordinate_string(p)
        handle.write(f'{chrom}\t{start}\t{end}\t{motif_length}\t{motif}\t{repeat_count}\n')
