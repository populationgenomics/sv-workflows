#!/usr/bin/env python3

"""
This script will output REViewer svg based on input (CPG ID) and locus ID, as defined in the variant catalog. 
Do not use this if inputting multiple sample IDs or locus IDs
"""
import os
import click


from cpg_utils.config import get_config
from cpg_utils.hail_batch import remote_tmpdir, output_path, reference_path
from cpg_workflows.batch import get_batch

REF_FASTA = 'gs://cpg-reference/hg38/v0/Homo_sapiens_assembly38.fasta'

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
REVIEWER_IMAGE = config['images']['reviewer']

@click.command()
@click.argument('id', )
@click.argument('locus')
@click.argument('catalog', 'GCP path to catalog')
@click.argument('input-dir', 'GCP path to input-dir, includes gs://')

def main(id:list[str], locus: str, catalog:str, input_dir: str):

    # Initializing Batch
    b= get_batch()
    samtools_job = b.new_job(name = f'Sorting and indexing')
    samtools_job.image(SAMTOOLS_IMAGE)
    samtools_job.storage('20G')
    samtools_job.cpu(8) 

    samtools_job.command(f"""

    """)

    var_inputs = b.read_input_group(**{'reads_bam': f'{input_dir}/'+id+ f'_EH_realigned_sorted.bam', 'vcf': f'{input_dir}/'+id +f'_EH.vcf',
                                       'reads_bai': f'{input_dir}/'+id+ f'_EH_realigned_sorted.bam.bai'})

    ref = b.read_input_group(
            **dict(
                base=REF_FASTA,
                catalog = catalog,
                fai=REF_FASTA + '.fai',
                dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
                + '.dict',
            )
        )  

    # Adding a job and giving it a descriptive name.
    reviewer_job = b.new_job(name = f'Visualising ' + locus)
    reviewer_job.image(REVIEWER_IMAGE)
    reviewer_job.storage('20G')
    reviewer_job.cpu(8)


    reviewer_job.declare_resource_group(ofile = {'svg': '{root}.'+locus+'.svg',
                                               'metrics.tsv': '{root}.metrics.tsv',
                                               'phasing.tsv': '{root}.phasing.tsv'
        })

    reviewer_job.command(f"""
        ./REViewer-v0.2.7-linux_x86_64 --reads {var_inputs['reads_bam']} --vcf {var_inputs['vcf']} --reference {ref.base} --catalog {ref.catalog} --out {reviewer_job.ofile} --locus {locus}
        """
        )
    # output writing    
    reviewer_out_fname = f'{locus}/{id}_{locus}_reviewer'
    reviewer_output_path = f'gs://cpg-tob-wgs-test/hoptan-str/tob10/reviewer/{reviewer_out_fname}'
    b.write_output(reviewer_job.ofile, reviewer_output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
    
    


