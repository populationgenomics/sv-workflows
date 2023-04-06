"""
utils script to prevent code repetiition
"""


from itertools import islice


def break_coordinate_string(coordinate_string: str) -> tuple[str, int, int]:
    """
    Break a coordinate string into chromosome, start, and end.
    :param coordinate_string:
    :return:
    """
    # split on the colon
    chrom, rest = coordinate_string.split(':')

    # split on the dash
    start, end = map(int, rest.split('-'))

    return chrom, start, end


def read_lines_from_file(file_handle, lines: int = 2):
    """
    Read N lines from a file handle, and return them as a list.
    """
    return list(islice(file_handle, lines))
