version 1.0

import "tasks/structs.wdl"

workflow Melt {
  input {
    File crbam
    File crbami
    String sample_id
    File fasta
    File? fastafai
    String? reference_version
    Float coverage
    Int read_length
    Float insert_size
    String docker_melt
    RuntimeAttr? runtime_attr_melt
  }

  call RunMELT {
    input:
      crbam = crbam,
      crbami = crbami,
      sample_id = sample_id,
      fasta = fasta,
      fastafai = fastafai,
      reference_version = reference_version,
      coverage = coverage,
      read_length = read_length,
      insert_size = insert_size,
      docker = docker_melt,
      runtime_attr_override = runtime_attr_melt
  }

  output {
    File vcf = RunMELT.vcf
    File vcfi = RunMELT.vcfi
  }
}

task RunMELT {
  input {
    File crbam
    File crbami
    String sample_id
    File fasta
    File? fastafai
    String? reference_version
    Float coverage
    Int read_length
    Float insert_size
    String docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    crbam: "Path to BAM/CRAM file (e.g. gs://path/to/sample.cram)."
    crbami: "Index for crbam (e.g. gs://path/to/sample.cram.crai)."
    sample_id: "Sample name, used for output (e.g. FOO will give a VCF output FOO.vcf.gz)"
    fasta: "Path to FASTA file used to align BAM or CRAM file (e.g. gs://ref/hg38.fasta)."
    fastafai: "[optional] The index for the FASTA file (default: '<reference_fasta>.fai') (e.g. gs://ref/hg38.fasta.fai)."
    reference_version: "[optional] Number of human genome assembly used (19 or 38) (default: 38)."
    coverage: "Median depth of coverage."
    read_length: "Read length."
    insert_size: "Insert size."
  }

  Boolean is_bam = basename(crbam, ".bam") + ".bam" == basename(crbam)
  File ref_index = select_first([fastafai, fasta + ".fai"])
  String ref_version = select_first([reference_version, "38"])

  # ensure there's sufficient disk space
  Float disk_overhead = 10.0
  Float crbam_size = size(crbam, "GiB")
  Float crbam_index_size = size(crbami, "GiB")
  Float ref_size = size(fasta, "GiB")
  Float ref_index_size = size(ref_index, "GiB")
  Int vm_disk_size = ceil(crbam_size + crbam_index_size + ref_size + ref_index_size + disk_overhead)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 16,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1,
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "${sample_id}.melt.vcf.gz"
    File vcfi = "${sample_id}.melt.vcf.gz.tbi"
  }

  command <<<
    set -Eeuo pipefail

    # MELT expects the BAM index to have extension ".bam.bai"
    mv ~{crbam} ~{sample_id}.cram
    mv ~{crbami} ~{sample_id}.cram.bai

    # these locations should be stable
    MELT_DIR="/MELT"
    CROMWELL_ROOT="/cromwell_root"

    # these locations may vary based on MELT version number, so find them:
    MELT_ROOT=$(find "${MELT_DIR}" -name "MELT.jar" | xargs -n1 dirname)
    MELT_SCRIPT="${MELT_DIR}/run_MELT.sh"

    # call MELT
    "${MELT_SCRIPT}" \
      "~{sample_id}.cram" \
      "~{fasta}" \
      ~{coverage} \
      ~{read_length} \
      ~{insert_size} \
      "$MELT_ROOT" \
      "$CROMWELL_ROOT" \
      ~{ref_version}

    # combine different mobile element VCFs into a single sample VCF
    # then sort into position order and compress with bgzip
    grep    "^#" SVA.final_comp.vcf   > "~{sample_id}.header.txt"
    grep -v "^#" SVA.final_comp.vcf   > "~{sample_id}.sva.vcf"
    grep -v "^#" LINE1.final_comp.vcf > "~{sample_id}.line1.vcf"
    grep -v "^#" ALU.final_comp.vcf   > "~{sample_id}.alu.vcf"
    cat ~{sample_id}.header.txt \
        ~{sample_id}.sva.vcf \
        ~{sample_id}.line1.vcf \
        ~{sample_id}.alu.vcf | \
      vcf-sort -c | \
      bgzip -c > ~{sample_id}.melt.vcf.gz
    tabix -p vcf "~{sample_id}.melt.vcf.gz"
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
