#!/usr/bin/env python3

"""

This script plots a histogram of the small DUPs genotyped by Wham in a single sample raw vcf.

analysis-runner --dataset "bioheart" --description "plot histogram of small DUPs" --access-level "test" \
    --output-dir "gatk_sv/qc_plots" \
    wham_small_dup_plotter.py --vcf-path gs://gatk-sv-ref-panel-1kg/outputs/tws-no-cram-conversion/GatherSampleEvidenceBatch/HG00096.wham.vcf.gz

"""

import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bokeh.plotting import output_file, save

import hail as hl

from cpg_utils.hail_batch import init_batch, output_path


@click.option('--vcf-path', help='GCS path to the input vcf file')
def main(vcf_path):
    init_batch()
    sample_id = vcf_path.split("/")[-1].split(".")[0]
    mt = hl.import_vcf(vcf_path, force=True)
    wham_tester_1_dup = mt.filter_rows(mt.info.SVTYPE == "DUP")
    wham_tester_1_dup = wham_tester_1_dup.annotate_rows(SVLEN_int=hl.int32(wham_tester_1_dup.info.SVLEN[0]))
    plot = hl.plot.histogram(
        wham_tester_1_dup.SVLEN_int,
        legend=f"Length of DUPs called by Wham in {sample_id} (raw vcf) (bp)",
        range=(0, 500),
    )
    output_file('local_plot.html')
    save(plot)
    gcs_path = output_path(f'wham_small_dups/{sample_id}.html', 'analysis')
    hl.hadoop_copy('local_plot.html', gcs_path)


if __name__ == '__main__':
    main()
