#!/usr/bin/env python3

"""
This script will output REViewer svg based on inputs: one/multiple CPG IDs and one locus, as defined in the variant catalog. 
analysis-runner --access-level test --dataset hgdp --description "reviewer" --output-dir 'str/410_sgd_loci/reviewer' reviewer_runner.py --catalog=gs://cpg-hgdp-test/str/410_sgdp_loci/catalogs/eh_catalog_hg38_backbone_trimmed_0_based.json --locus=chr3-67969584-67969611-AAC --input-dir=gs://cpg-hgdp-test/str/sensitivity-analysis/eh/trimmed_coordinates_0_based_hg38_backbone CPG265538		

"""
import os
import click


from cpg_utils.config import get_config
from cpg_utils.hail_batch import output_path
from cpg_workflows.batch import get_batch

REF_FASTA = 'gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta'

config = get_config()

SAMTOOLS_IMAGE = config['images']['samtools']
REVIEWER_IMAGE = config['images']['reviewer']

@click.option('--locus', help = 'Locus identifier as per catalog')
@click.option('--catalog', help='GCP path to catalog')
@click.option('--input-dir', help='GCP path to input-dir, includes gs://')
@click.argument('cpg-wgs-ids',nargs =-1 )
@click.command()


def main(locus,catalog,input_dir,cpg_wgs_ids:list[str]):

    # Initializing Batch
    b= get_batch()
    ref = b.read_input_group(
                **dict(
                    base=REF_FASTA,
                    catalog = catalog,
                    fai=REF_FASTA + '.fai',
                    dict=REF_FASTA.replace('.fasta', '').replace('.fna', '').replace('.fa', '')
                    + '.dict',
                )
            )  
    for cpg_id in cpg_wgs_ids: 
        #read in bam file from EH output
        bam_input = b.read_input( f'{input_dir}/{cpg_id}_eh.realigned_bam')

        #read in VCF from EH output
        vcf_input = b.read_input( f'{input_dir}/{cpg_id}_eh.vcf')

        #BAM file must be sorted and indexed prior to inputting into REViewer 
        samtools_job = b.new_job(name = f'Sorting and indexing {cpg_id}')
        samtools_job.image(SAMTOOLS_IMAGE)
        samtools_job.storage('20G')
        samtools_job.cpu(8) 
        samtools_job.declare_resource_group(
            bam= {
            'sorted.bam':'{root}.sorted.bam',
            'sorted.bam.bai':'{root}.sorted.bam.bai'
            }
        )

        samtools_job.command(f"""

        echo "sorting {bam_input}";
        samtools sort {bam_input} -o {samtools_job.bam['sorted.bam']};

        echo "indexing {samtools_job.bam['sorted.bam']}";
        samtools index {samtools_job.bam['sorted.bam']} -b -o {samtools_job.bam['sorted.bam.bai']};

        """)

        # REViewer job 
        reviewer_job = b.new_job(name = f'Visualising {locus} for {cpg_id}')
        reviewer_job.image(REVIEWER_IMAGE)
        reviewer_job.depends_on(samtools_job)
        reviewer_job.storage('20G')
        reviewer_job.cpu(8)
        
        reviewer_job.declare_resource_group(ofile = {'svg': '{root}.'+f'{locus}.svg',
                                               'metrics.tsv': '{root}.metrics.tsv',
                                               'phasing.tsv': '{root}.phasing.tsv'
        })

        reviewer_job.command(f"""
            ./REViewer-v0.2.7-linux_x86_64 --reads {samtools_job.bam['sorted.bam']} --vcf {vcf_input} --reference {ref.base} --catalog {ref.catalog} --out {reviewer_job.ofile} --locus {locus}
            """
            )
    # output writing    
    reviewer_output_path = output_path(f'{locus}/{cpg_id}_{locus}_reviewer', 'analysis')
    b.write_output(reviewer_job.ofile, reviewer_output_path)

    b.run(wait=False)

if __name__ == '__main__':
    main() 
    
    


