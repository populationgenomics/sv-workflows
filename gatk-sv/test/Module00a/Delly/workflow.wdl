version 1.0

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.14-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.14-beta/wdl/Delly.wdl" as delly

workflow Delly {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File exclude_intervals_file
    String sv_base_mini_docker
    String delly_docker
  }
  Array[String] event_types = ["DEL", "DUP", "INV"]

  scatter (event_type in event_types) {
    call delly.RunDelly {
      input:
        bam_or_cram_file = bam_or_cram_file,
        bam_or_cram_index = bam_or_cram_index,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        exclude_intervals_file = exclude_intervals_file,
        event_type = event_type,
        delly_docker = delly_docker,
    }
  }

  call delly.GatherBCFs {
    input:
      bcfs = RunDelly.bcf,
      sample_id = sample_id,
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File vcf = GatherBCFs.vcf
    File index = GatherBCFs.index
  }
}
