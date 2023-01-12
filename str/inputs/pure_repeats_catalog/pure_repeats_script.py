"""
script alternative to the ipython notebook code
"""


TRF_DATA = 'trf_output.dat'
BED_CATALOG = 'bed_catalog_without_complex_repeats.bed'
SEQUENCE = 'Sequence: '
PARAMETERS = 'Parameters:'
FASTA_START = '>'

# indices
PURITY_INDEX = 5
REPEAT_INDEX = 13


def get_purest_motif(metric_list: list[str]) -> tuple[str, int]:
    """
    take the list of metrics lines, find the purest
    wasteful if there's only one row
    allows for multiple motifs at one locus
    :param metric_list:
    :return: the purest motif and assc. purity
    """
    metric_data = []
    for each_metric in metric_list:
        split_line = each_metric.rstrip().split()
        metric_data.append((split_line[REPEAT_INDEX], int(split_line[PURITY_INDEX])))

    return max(metric_data, key=lambda x: x[1])


sections = []

# parse the file into sections, broken up by 'Sequence:' lines
# each 'section' will contain the sequence, params, and metrics lines
with open(TRF_DATA, encoding='utf-8') as filehandle:

    # start a fresh section
    this_section: list[str] = []

    # iterate over lines in this file
    for line in filehandle:

        # demarcates a new sequence, so add the old
        if line.startswith(SEQUENCE):
            # avoid storing the header
            if any(section_line.startswith(SEQUENCE) for section_line in this_section):
                sections.append(this_section)

            # reset the section - adding this line
            this_section = [line.rstrip()]

        # ignore blank lines, but keep adding to the existing section
        elif line.rstrip():
            this_section.append(line.rstrip())

    # after iterating over the file, catch the final section
    sections.append(this_section)


blank_coords = []
pure_repeats = []

# process the sections
for section in sections:

    # if there were no metrics - get the blank coordinates
    if len(section) == 2:
        blank_coords.append(section[0].removeprefix(SEQUENCE).rstrip())
        continue

    # otherwise get the purest repeat, and the purity score
    # not interested in parameters?
    this_seq = section[0].removeprefix(SEQUENCE).rstrip()
    motif, purity = get_purest_motif(section[2:])
    if purity == 100:
        pure_repeats.append((this_seq, motif))


# debug
print(len(blank_coords))  # 6409
print(len(pure_repeats))  # 161612


# bed catalog bits
motif_dictionary = {}
with open('bed_catalog_without_complex_repeats.bed', encoding='utf-8') as handle:
    for line in handle:
        line_list = line.split('\t')
        motif_dictionary[f'{line_list[0]}:{line_list[1]}-{line_list[2]}'] = line_list[
            4
        ].rstrip()


def chunks(iterable, chunk_size):
    """
    Yield successive n-sized chunks from an iterable
    """

    for i in range(0, len(iterable), chunk_size):
        yield iterable[i : (i + chunk_size)]


# fasta_bits
with open('catalog_fasta_sequences.fasta.txt', encoding='utf-8') as handle:
    fasta_lines = handle.readlines()

fasta_dict = {}
# iterate over lines in pairs. not efficient, but whatever
# would need re-thinking if the motifs stretched to multiple lines
for line_pair in chunks(fasta_lines, 2):
    assert line_pair[0].startswith(FASTA_START)
    fasta_dict[line_pair[0].lstrip(FASTA_START).rstrip()] = line_pair[1].rstrip()


# debug
print(fasta_dict['chr12:6936716-6936773'])
print(motif_dictionary['ßchr12:6936716-6936773'])

# check if any of the impure repeats can be recovered
impure_repeats = []
for blank in blank_coords:
    motif = motif_dictionary[blank]
    fasta = fasta_dict[blank]
    num_repeats = len(fasta) / len(motif)
    if fasta != motif * int(num_repeats):
        impure_repeats.append(blank)
    else:
        pure_repeats.append((blank, motif))

print(len(pure_repeats))

with open('output.txt', 'w', encoding='utf-8') as handle:
    handle.writelines('\n'.join(['\t'.join(pure) for pure in pure_repeats]))
