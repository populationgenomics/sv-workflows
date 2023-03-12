#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,no-member

import os
import math

from cloudpathlib import AnyPath
from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch

REF_FASTA = 'gs://hg19_hipstr_experiment/hs37d5.fa.gz'

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
REVIEWER_IMAGE = config['images']['reviewer']

def main(): 
    b = get_batch()
    ref = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
    )
    samtools_job = b.new_job("Indexer")
    samtools_job.image(SAMTOOLS_IMAGE)

    samtools_job.declare_resource_group(
            bam={
                'sorted.bam': '{root}.sorted.bam',
                'sorted.bam.bai': '{root}.sorted.bam.bai',
            }
        )
    bam_input = b.read_input("gs://hg19_hipstr_experiment/LP6005441-DNA_A01.bam")
    samtools_job.cpu(4)  # to support the -@ 4
    samtools_job.command(
            f"""
        echo "sorting {bam_input}";
        samtools sort -@ 4 {bam_input} -o {samtools_job.bam['sorted.bam']};
        echo "indexing {samtools_job.bam['sorted.bam']}";
        samtools index -@ 4 {samtools_job.bam['sorted.bam']} -b -o {samtools_job.bam['sorted.bam.bai']};
        """
        )
    # output writing
    samtools_output_path = output_path(
            'samtools'
        )
    b.write_output(samtools_job.bam, samtools_output_path)

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,missing-function-docstring
