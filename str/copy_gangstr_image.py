#!/usr/bin/env python


import os
import hailtop.batch as hb
import click

DATASET = os.getenv('DATASET')
HAIL_BUCKET = os.getenv('HAIL_BUCKET')
OUTPUT_SUFFIX = os.getenv('OUTPUT')
BILLING_PROJECT = os.getenv('HAIL_BILLING_PROJECT')
ACCESS_LEVEL = os.getenv('ACCESS_LEVEL')

REF_FASTA = 'gs://cpg-reference/hg38/v1/Homo_sapiens_assembly38.fasta'
SAMTOOLS_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/samtools:v0'


@click.command()

def main():  # pylint: disable=missing-function-docstring
    """
    Subset CRAM or BAM file CRAM_PATH to REGION. Example: batch.py sample.cram chr21:1-10000
    """
    # Initializing Batch
    backend = hb.ServiceBackend(billing_project=BILLING_PROJECT, bucket=HAIL_BUCKET)
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    # Adding a job and giving it a descriptive name.
    j = b.new_job('Pulling GangSTR image')


    # This image contains basic bioinformatics tools like samtools, bcftools, Picard, etc.
    j.image(SAMTOOLS_IMAGE)

    # For larger CRAMs, request more storage.
    j.storage('10G')

    # If you want to run a multithreaded command, e.g. samtools with -@, request more CPUs here.
    # Note that the machines use hyperthreading, so for every CPU, 2x threads are available.
    j.cpu(2)

    # The command that do the actual job.
    j.command(
        f"""
        set -ex

        gcloud auth configure-docker australia-southeast1-docker.pkg.dev
        docker pull weisburd/gangstr:v2.5
        docker tag weisburd/gangstr:v2.5 australia-southeast1-docker.pkg.dev/cpg-common/images/gangstr:v2.5
        docker push australia-southeast1-docker.pkg.dev/cpg-common/images/gangstr:v2.5   
    """
    )


#skopeo copy docker://weisburd/gangstr:v2.5 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/gangstr:v2.5


    # don't wait for the hail batch workflow to complete, otherwise
    # the workflow might get resubmitted if this VM gets preempted.
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter