version 1.0

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/CollectCoverage.wdl" as cov

workflow CollectCoverage {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File preprocessed_intervals
    File reference_fasta
    File reference_index
    File reference_dict
    String gatk_docker
  }

  call cov.CollectCounts {
    input:
      bam = bam_or_cram_file,
      bam_idx = bam_or_cram_index,
      sample_id = sample_id,
      intervals = preprocessed_intervals,
      ref_fasta = reference_fasta,
      ref_fasta_fai = reference_index,
      ref_fasta_dict = reference_dict,
      gatk_docker = gatk_docker,
      disabled_read_filters = ["MappingQualityReadFilter"]
  }
}
