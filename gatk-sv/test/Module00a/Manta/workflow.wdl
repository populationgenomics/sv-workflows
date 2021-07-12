version 1.0

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Manta.wdl" as manta

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
