#!/usr/bin/env python3

"""
Calculates GC metrics for a CRAM file using Picard CollectGcBiasMetrics.

analysis-runner --dataset "bioheart" \
    --description "Compute GC metrics" \
    --access-level "test" \
    --output-dir "qc-stand-alone" \
    compute_gc_metrics.py --input-dir "gs://cpg-bioheart-test/cram" CPGXXX CPGYYY


"""
import os
import click

from cpg_utils.hail_batch import command, image_path, get_batch, output_path
from cpg_workflows.resources import HIGHMEM


REF_FASTA = 'gs://cpg-common-main/references/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta'


# input directory to CRAMs
@click.option('--input-dir', help='gs://...')
# input CPG sample ID
@click.argument('internal-wgs-ids', nargs=-1)
@click.command()
def main(input_dir, internal_wgs_ids):
    """
    Make job that runs Picard CollectGcBiasMetrics.
    """
    b = get_batch(name='Compute GC metrics')
    reference = b.read_input_group(
        base=REF_FASTA,
        fai=REF_FASTA + '.fai',
        dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
        + '.dict',
    )

    cram_dict = {id: os.path.join(input_dir, f'{id}.cram') for id in internal_wgs_ids}

    for id in list(cram_dict.keys()):

        input_cram = b.read_input(cram_dict[id])

        j = b.new_job(f'CollectGcBiasMetrics: {id}')

        j.image(image_path('picard'))
        resource = HIGHMEM.request_resources(ncpu=4)
        resource.attach_disk_storage_gb = 250
        resource.set_to_job(j)

        j.declare_resource_group(
            output_files={
                'gc_bias_sample.txt': '{root}.gc_bias_metrics.txt',
                'gc_bias_summary.txt': '{root}.summary_metrics.txt',
                'gc_bias_chart.pdf': '{root}.gc_bias_metrics.pdf',
            }
        )

        cmd = f"""


        picard CollectGcBiasMetrics -Xms{resource.get_java_mem_mb()}M \\
        I={input_cram} O={j.output_files['gc_bias_sample.txt']} CHART={j.output_files['gc_bias_chart.pdf']} \\
        S={j.output_files['gc_bias_summary.txt']} \\
        R={reference.base}
        echo "CollectGcBiasMetrics finished successfully"

        """
        j.command(command(cmd, monitor_space=True))

        b.write_output(j.output_files, output_path(f'{id}', 'analysis'))

    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter,missing-function-docstring
