version 1.0

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/CollectCoverage.wdl" as cov
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/CramToBam.wdl" as ctb
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Delly.wdl" as delly
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Manta.wdl" as manta
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/PESRCollection.wdl" as pesr
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Whamg.wdl" as wham

workflow Module00a {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index
    String sample_id

    # Use only for crams in requester pays buckets
    Boolean requester_pays_crams = false

    # Evidence collection flags
    Boolean collect_coverage = true
    Boolean collect_pesr = true

    # Common parameters
    File primary_contigs_list
    File reference_fasta
    File reference_index    # Index (.fai), must be in same dir as fasta
    File reference_dict     # Dictionary (.dict), must be in same dir as fasta

    # Coverage collection inputs
    File preprocessed_intervals

    # Delly inputs
    File? delly_exclude_intervals_file  # Required if run_delly True

    # Manta inputs
    File manta_region_bed
    File? manta_region_bed_index

    # Wham inputs
    File wham_include_list_bed_file

    # Docker
    String sv_base_mini_docker
    String samtools_cloud_docker
    String? delly_docker
    String? manta_docker
    String? wham_docker
    String gatk_docker
    String? gatk_docker_pesr_override

    # Never assign these values! (workaround until None type is implemented)
    Float? NONE_FLOAT_
    Int? NONE_INT_
    File? NONE_FILE_
  }

  Boolean run_delly = defined(delly_docker)
  Boolean run_manta = defined(manta_docker)
  Boolean run_wham = defined(wham_docker)

  Boolean is_bam_ = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  String index_ext_ = if is_bam_ then ".bai" else ".crai"
  File bam_or_cram_index_ = if defined(bam_or_cram_index) then select_first([bam_or_cram_index]) else bam_or_cram_file + index_ext_

  # Convert to BAM if we have a CRAM
  if (!is_bam_) {
    call ctb.CramToBam {
      input:
        cram_file = bam_or_cram_file,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        requester_pays = requester_pays_crams,
        samtools_cloud_docker = samtools_cloud_docker,
    }
  }

  File bam_file_ = select_first([CramToBam.bam_file, bam_or_cram_file])
  File bam_index_ = select_first([CramToBam.bam_index, bam_or_cram_index_])

  if (collect_coverage) {
    call cov.CollectCounts {
      input:
        intervals = preprocessed_intervals,
        bam = bam_file_,
        bam_idx = bam_index_,
        sample_id = sample_id,
        ref_fasta = reference_fasta,
        ref_fasta_fai = reference_index,
        ref_fasta_dict = reference_dict,
        gatk_docker = gatk_docker,
        disabled_read_filters = ["MappingQualityReadFilter"]
    }
  }

  if (run_delly) {
    call delly.Delly {
      input:
        bam_or_cram_file = bam_file_,
        bam_or_cram_index = bam_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        exclude_intervals_file = select_first([delly_exclude_intervals_file]),
        sv_base_mini_docker = sv_base_mini_docker,
        delly_docker = select_first([delly_docker]),
    }
  }

  if (run_manta) {
    call manta.Manta {
      input:
        bam_or_cram_file = bam_file_,
        bam_or_cram_index = bam_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        region_bed = manta_region_bed,
        manta_docker = select_first([manta_docker]),
    }
  }

  if (collect_pesr) {
    call pesr.PESRCollection {
      input:
        cram = bam_file_,
        cram_index = bam_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        reference_dict = reference_dict,
        gatk_docker = select_first([gatk_docker_pesr_override, gatk_docker]),
    }
  }

  if (run_wham) {
    call wham.Whamg {
      input:
        bam_or_cram_file = bam_file_,
        bam_or_cram_index = bam_index_,
        sample_id = sample_id,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        include_bed_file = wham_include_list_bed_file,
        chr_file = primary_contigs_list,
        samtools_cloud_docker = samtools_cloud_docker,
        wham_docker = select_first([wham_docker]),
    }
  }

  output {
    File? coverage_counts = CollectCounts.counts

    File? delly_vcf = Delly.vcf
    File? delly_index = Delly.index

    File? manta_vcf = Manta.vcf
    File? manta_index = Manta.index

    File? pesr_disc = PESRCollection.disc_out
    File? pesr_disc_index = PESRCollection.disc_out_index
    File? pesr_split = PESRCollection.split_out
    File? pesr_split_index = PESRCollection.split_out_index

    File? wham_vcf = Whamg.vcf
    File? wham_index = Whamg.index
  }
}
