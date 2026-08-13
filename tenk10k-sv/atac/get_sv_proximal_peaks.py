#!/usr/bin/env python3
# pylint: disable=too-many-locals

"""
This script aims to:
 - identify the ATAC-seq peaks ('loc' column of PseudobulkMatStandardised.csv) that fall within
   a window (default 1Mb) of any SV in a VCF, per cell type.
 - output {cell_type}_loc.csv (a single 'loc' column), which is the format expected by the
   loc_filtered_dir input of str/atac/get_cis_numpy.py.

SV coordinates are taken as POS to INFO/END; SVs without an END (eg BND, INS) are treated as
POS to POS. A peak is retained if it overlaps [POS - window_size, END + window_size].

analysis-runner --config get_sv_proximal_peaks.toml --dataset "bioheart" --access-level "test" \
--description "get SV-proximal ATAC peaks" --output-dir "str/atac-seq/sv-proximal-peaks/1mb" \
python3 get_sv_proximal_peaks.py

"""

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, image_path, output_path


def sv_proximal_loc_extractor(
    sv_intervals_file,
    pseudobulk_dir,
    cell_type,
    window_size,
):
    """
    Writes the subset of peaks lying within window_size of an SV to {cell_type}_loc.csv

    """
    import numpy as np
    import pandas as pd

    from cpg_utils.hail_batch import output_path

    def strip_chr(series):
        # the pseudobulk matrices use 'chr1', so normalise both sides rather than assume the VCF matches
        return series.astype(str).str.replace('^chr', '', regex=True)

    svs = pd.read_csv(sv_intervals_file, sep='\t', names=['chrom', 'pos', 'end'], dtype=str)
    svs['pos'] = svs['pos'].astype(int)
    # bcftools writes '.' when INFO/END is absent
    svs['end'] = pd.to_numeric(svs['end'], errors='coerce').fillna(svs['pos']).astype(int)
    svs['chrom'] = strip_chr(svs['chrom'])
    print(f'Read {svs.shape[0]} SVs across {svs["chrom"].nunique()} contigs')

    peaks = pd.read_csv(f'{pseudobulk_dir}/{cell_type}/PseudobulkMatStandardised.csv', usecols=['loc'])
    coords = peaks['loc'].astype(str).str.extract(r'^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$')
    unparsed = coords['chrom'].isna().sum()
    if unparsed:
        print(f'WARNING: dropping {unparsed} peaks whose loc could not be parsed as chrom:start-end')
    peaks = peaks.assign(
        chrom=strip_chr(coords['chrom']),
        start=pd.to_numeric(coords['start']),
        end=pd.to_numeric(coords['end']),
    ).dropna(subset=['start', 'end'])
    print(f'Read {peaks.shape[0]} peaks for {cell_type}')

    shared_contigs = set(peaks['chrom']) & set(svs['chrom'])
    if not shared_contigs:
        raise ValueError(
            f'No contigs shared between the SV VCF and the {cell_type} peaks - check chromosome naming. '
            f'SV contigs: {sorted(set(svs["chrom"]))[:5]}, peak contigs: {sorted(set(peaks["chrom"]))[:5]}',
        )

    retained = []
    for chrom, chrom_peaks in peaks.groupby('chrom'):
        chrom_svs = svs[svs['chrom'] == chrom]
        if chrom_svs.empty:
            continue

        window_starts = np.maximum(chrom_svs['pos'].to_numpy() - window_size, 0)
        window_ends = chrom_svs['end'].to_numpy() + window_size
        order = np.argsort(window_starts)
        window_starts, window_ends = window_starts[order], window_ends[order]

        # merge the windows into disjoint intervals so each peak can be tested with a single lookup
        merged_starts: list[int] = []
        merged_ends: list[int] = []
        for start, end in zip(window_starts, window_ends):
            if merged_ends and start <= merged_ends[-1]:
                merged_ends[-1] = max(merged_ends[-1], end)
            else:
                merged_starts.append(start)
                merged_ends.append(end)
        merged_starts_arr = np.array(merged_starts)
        merged_ends_arr = np.array(merged_ends)

        peak_starts = chrom_peaks['start'].to_numpy()
        peak_ends = chrom_peaks['end'].to_numpy()
        # the last window starting at or before the peak end is the only one that can overlap
        idx = np.searchsorted(merged_starts_arr, peak_ends, side='right') - 1
        overlaps = (idx >= 0) & (merged_ends_arr[np.clip(idx, 0, None)] >= peak_starts)
        retained.append(chrom_peaks[overlaps])

    proximal = pd.concat(retained) if retained else peaks.iloc[:0]
    print(
        f'{cell_type}: {proximal.shape[0]} of {peaks.shape[0]} peaks '
        f'({100 * proximal.shape[0] / max(peaks.shape[0], 1):.1f}%) within {window_size}bp of an SV',
    )

    proximal[['loc']].to_csv(output_path(f'{cell_type}_loc.csv'), index=False)


def main():
    """
    Extract SV intervals once, then identify SV-proximal peaks per cell type
    """
    b = get_batch(name='get SV-proximal ATAC peaks')
    config = get_config()['get_sv_proximal_peaks']

    # pull CHROM/POS/END out of the VCF once, so the per-cell-type jobs only handle a small table
    extract_job = b.new_job(name='Extract SV intervals')
    extract_job.image(image_path('bcftools'))
    extract_job.storage(config['extract_job_storage'])
    vcf = b.read_input(config['sv_vcf_path'])
    extract_job.command(
        f"bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\n' {vcf} > {extract_job.ofile}",
    )
    b.write_output(extract_job.ofile, output_path('sv_intervals.tsv'))

    for cell_type in config['cell_types'].split(','):
        j = b.new_python_job(name=f'Identify SV-proximal peaks: {cell_type}')
        j.cpu(config['job_cpu'])
        j.memory(config['job_memory'])
        j.storage(config['job_storage'])

        j.call(
            sv_proximal_loc_extractor,
            extract_job.ofile,
            config['pseudobulk_dir'],
            cell_type,
            config['window_size'],
        )

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
