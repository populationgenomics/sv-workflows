#!/usr/bin/env python3

"""
Running Batch script to produce QC output for GangSTR calls using TRTools. 
"""
import os
import hailtop.batch as hb
import click

DATASET = os.getenv('DATASET')
HAIL_BUCKET = os.getenv('HAIL_BUCKET')
OUTPUT_SUFFIX = os.getenv('OUTPUT')
BILLING_PROJECT = os.getenv('HAIL_BILLING_PROJECT')
ACCESS_LEVEL = os.getenv('ACCESS_LEVEL')


TRTOOLS_IMAGE = "australia-southeast1-docker.pkg.dev/cpg-common/images/trtools:v4.0.2"
BCFTOOLS_IMAGE = "australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools:1.10.2--h4f4756c_2"
expansionhunter_vcf_path = "gs://cpg-tob-wgs-test/hoptan-str/tob10/merge/tob10_EH_merge.vcf"

@click.command()
def main():  # pylint: disable=missing-function-docstring

    # Initializing Batch
    backend = hb.ServiceBackend(billing_project=BILLING_PROJECT, bucket=HAIL_BUCKET)
    b = hb.Batch(backend=backend, default_image=os.getenv('DRIVER_IMAGE'))

    vcf = b.read_input_group(**{'eh': expansionhunter_vcf_path})

    # BCFTools is needed to sort, zip, and index files before using them as input in TRTools 
    bcftools_job = b.new_job("Files prep")
    bcftools_job.image(BCFTOOLS_IMAGE)
    bcftools_job.storage('200G')
    bcftools_job.cpu(8)

    # declare resource groups, including extensions
    bcftools_job.declare_resource_group(
        eh_vcf={'vcf.gz': '{root}.eh.vcf.gz', 'vcf.gz.tbi': '{root}.eh.vcf.gz.tbi'}
    )

    bcftools_job.command(f"""
    set -ex;
    
    echo "compressing {vcf['eh']}"; 
    bgzip -c {vcf['eh']} > {bcftools_job.eh_vcf['vcf.gz']};
    
    echo "indexing {bcftools_job.eh_vcf['vcf.gz']}";
    tabix -p vcf {bcftools_job.eh_vcf['vcf.gz']};
    """)

    #TRTools package contains qcSTR() to output QC graphs. 
    trtools_job = b.new_job("Run qcSTR")
    trtools_job.image(TRTOOLS_IMAGE)
    trtools_job.depends_on(bcftools_job)

    trtools_job.storage('200G')
    trtools_job.cpu(8)

    trtools_job.declare_resource_group(ofile = {
                                               'diffref-bias-modified.pdf': '{root}-diffref-bias.pdf'
                                                #'quality.pdf': '{root}-quality.pdf' EH does not have quality plots 

    })

    trtools_job.command(f"""
    set -ex;
    qcSTR --vcf {bcftools_job.eh_vcf['vcf.gz']} --vcftype eh  --refbias-binsize 1 --refbias-xrange-min 39 --out {trtools_job.ofile}
    
    """
    )
    # Speciying where to write the result
    out_fname = 'qcSTR_EH'
    output_path = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/qcSTR/eh/{out_fname}'
    b.write_output(trtools_job.ofile, output_path)

    # don't wait for the hail batch workflow to complete, otherwise
    # the workflow might get resubmitted if this VM gets preempted.
    b.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
    