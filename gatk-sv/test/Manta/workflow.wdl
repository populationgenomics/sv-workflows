version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.15-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.15-beta/wdl/Manta.wdl" as manta

workflow Manta {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File region_bed
    File? region_bed_index
    Float? jobs_per_cpu
    Int? mem_gb_per_job
    String manta_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk."
    bam_or_cram_index: "[optional] Index for bam_or_cram_file. If not specified, index is assumed to be at bam_file_path + '.bai' or cram_file_path + '.crai'"
    sample_id: "sample name. Outputs will be sample_name + 'manta.vcf.gz' and sample_name + 'manta.vcf.gz.tbi'"
    reference_fasta: ".fasta file with reference used to align bam or cram file"
    reference_index: "[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file"
    region_bed: "[optional] gzipped bed file with included regions where manta should make SV calls."
    region_bed_index: "[optional]If omitted, the WDL will look for an index by appending .tbi to the region_bed file"
    jobs_per_cpu: "[optional] number of manta threads, i.e. num_jobs = round(num_cpu * jobs_per_cpu). If omitted, defaults to 1.3."
    mem_gb_per_job: "[optional] Memory to request for VM (in GB) = round(num_jobs * mem_gb_per_job). If omitted, defaults to 2."
  }

  call manta.RunManta {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      sample_id = sample_id,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      region_bed = region_bed,
      region_bed_index = region_bed_index,
      jobs_per_cpu = jobs_per_cpu,
      mem_gb_per_job = mem_gb_per_job,
      manta_docker = manta_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File vcf = RunManta.vcf
    File index = RunManta.index
  }
}
