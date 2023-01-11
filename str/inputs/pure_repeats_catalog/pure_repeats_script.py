"""
script alternative to the ipython notebook code
"""


import sys
import csv
import json
import re


TRF_DATA = 'trf_output.dat'
SEQUENCE = 'Sequence:'
PARAMETERS = 'Parameters:'

# indices
PURITY_INDEX = 5
REPEAT_INDEX = 13


def get_purest_motif(metric_list: list[str]) -> tuple[str, int]:
    """
    take the list of metrics lines, find the purest
    wasteful if there's only one row
    :param metric_list:
    :return:
    """
    metric_data = []
    for each_metric in metric_list:
        line_list = each_metric.rstrip().split()
        metric_data.append((line_list[REPEAT_INDEX], int(line_list[PURITY_INDEX])))

    return max(metric_data, key=lambda x: x[1])


sections = []

# parse the file into sections, broken up by 'Sequence:' lines
with open(sys.argv[1]) as filehandle:

    this_section = []
    for line in filehandle:
        if line.startswith(SEQUENCE):
            # avoid storing the header
            if any(section_line.startswith(SEQUENCE) for section_line in this_section):
                sections.append(this_section)

            # reset the section
            this_section = []

        # ignore blank lines
        if line.rstrip():
            this_section.append(line.rstrip())

    sections.append(this_section)

blank_coords = []
pure_repeats = []

# process the sections:
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
        pure_repeats.append((this_seq, motif, purity))


print(len(blank_coords))  # 6409
print(len(pure_repeats))  # 161612

