version 1.0

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.14-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.14-beta/wdl/PESRCollection.wdl" as PESRC

# Workflow to run PE/SR collection on a single sample
workflow PESRCollection {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File reference_index
    File reference_dict
    String gatk_docker
  }

  call PESRC.RunPESRCollection {
    input:
      cram = bam_or_cram_file,
      cram_index = bam_or_cram_index,
      sample_id = sample_id,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      reference_dict = reference_dict,
      gatk_docker = gatk_docker
  }

  output {
    File disc_out = RunPESRCollection.disc_out
    File disc_out_index = RunPESRCollection.disc_out_index
    File split_out = RunPESRCollection.split_out
    File split_out_index = RunPESRCollection.split_out_index
  }
}
