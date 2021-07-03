version 1.0

import "../../gatk-sv-git/wdl/Whamg.wdl" as wham
import "../../gatk-sv-git/wdl/Structs.wdl"
import "../../gatk-sv-git/wdl/CramToBam.wdl" as ctb

# Run Whamg SV detection algorithm on whole genome in bam or cram
#   file, or if include_list is provided, run whamg on explicitly included
#   subset of genome and concatenate.
# Then fix output vcf headers to contain good contig length data.

workflow Whamg {
  input {
    File bam_or_cram_file
    File? bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File? include_bed_file
    File chr_file
    Int? pf_reads_improper_pairs
    Float? pct_exc_total
    String samtools_cloud_docker
    String wham_docker
    RuntimeAttr? runtime_attr_includelist
    RuntimeAttr? runtime_attr_cram_to_bam
    RuntimeAttr? runtime_attr_wham
  }
    
  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. bams are preferable, crams will be converted to bams."
    bam_or_cram_index: "[optional] associated index file. If omitted, the WDL will look for an index file by appending .bai/.crai to the .bam/.cram file"
    sample_id: "sample name. Outputs will be sample_name + '.vcf.gz' and sample_name + '.vcf.gz.tbi'"
    reference_fasta: ".fasta file with reference used to align bam or cram file"
    reference_index: "[optional] reference index file. If omitted, the WDL will look for an index by appending .fai to the .fasta file"
    chr_file: "text file with newline-separated list of contigs that whamg will use to estimate template size. Typically you will want only primary contigs."
    include_bed_file: "[optional] bed file with intervals where whamg should make calls. If omitted, whamg will run on whole genome."
    pf_reads_improper_pairs: "[optional] Value of PF_READS_IMPROPER_PAIRS obtained from CollectAlignmentSummaryMetrics, used for optimal estimate of VM memory needs."
    pct_exc_total: "[optional] Value of PCT_EXC_TOTAL obtained from CollectWgsMetrics, used for optimal estimate of VM memory needs."
  }

  if (basename(bam_or_cram_file, ".bam") + ".bam" != basename(bam_or_cram_file)) {
    call ctb.RunCramToBam as RunCramToBam {
      input:
        cram_file = bam_or_cram_file,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        samtools_cloud_docker = samtools_cloud_docker,
        runtime_attr_override = runtime_attr_cram_to_bam
    }
  }

  File bam_file = select_first([RunCramToBam.bam_file, bam_or_cram_file])
  File bam_index = select_first([RunCramToBam.bam_index, bam_or_cram_index, bam_or_cram_file + ".bai"])
  
  # decide whether to use includelist version or baseline
  Boolean use_include_list = defined(include_bed_file)

  if (use_include_list) {
    call wham.RunWhamgIncludelist {
      input:
        bam_file = bam_file,
        bam_index = bam_index,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        sample_id = sample_id,
        chr_file = chr_file,
        include_bed_file = select_first([include_bed_file]),
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        pct_exc_total = pct_exc_total,
        wham_docker = wham_docker,
        runtime_attr_override = runtime_attr_includelist
    }
  } # else
  if (!use_include_list) {
    call wham.RunWhamg {
      input:
        bam_file = bam_file,
        bam_index = bam_index,
        reference_fasta = reference_fasta,
        reference_index = reference_index,
        sample_id = sample_id,
        chr_file = chr_file,
        pf_reads_improper_pairs = pf_reads_improper_pairs,
        pct_exc_total = pct_exc_total,
        wham_docker = wham_docker,
        runtime_attr_override = runtime_attr_wham
    }
  }

  output {
    File index = select_first([RunWhamg.index, RunWhamgIncludelist.index])
    File vcf = select_first([RunWhamg.vcf, RunWhamgIncludelist.vcf])
  }
}
