version 1.0

import "structs.wdl"

workflow Delly {
  input {
    File crbam
    File? crbami
    String sample_id
    File fasta
    File? fastafai
    File blocklist
    Array[String]? sv_types
    String docker_sv_mini_base
    String docker_delly
    RuntimeAttr? runtime_attr_delly
    RuntimeAttr? runtime_attr_gather
  }

  Array[String] sv_types_array = select_first([sv_types, ["DEL", "DUP", "INV"]])

  scatter (sv_type in sv_types_array) {
    call RunDelly {
      input:
        crbam = crbam,
        crbami = crbami,
        sample_id = sample_id,
        fasta = fasta,
        fastafai = fastafai,
        blocklist = blocklist,
        sv_type = sv_type,
        docker = docker_delly,
        runtime_attr_override = runtime_attr_delly
    }
  }

  call GatherBCFs {
    input:
      bcfs = RunDelly.bcf,
      sample_id = sample_id,
      docker = docker_sv_mini_base,
      runtime_attr_override = runtime_attr_gather
  }

  output {
    File vcf = GatherBCFs.vcf
    File index = GatherBCFs.index
  }
}

task RunDelly {
  input {
    File crbam
    File? crbami
    String sample_id
    File fasta
    File? fastafai
    File blocklist
    String sv_type
    String docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    crbam: "Path to BAM/CRAM file (e.g. gs://path/to/sample.cram)."
    crbami: "[optional] Index for crbam (default: '<crbam>.bai/.crai') (e.g. gs://path/to/sample.cram.crai)."
    sample_id: "Sample name, used for output (e.g. FOO will give a VCF output FOO.delly.vcf.gz)"
    fasta: "Path to FASTA file used to align BAM or CRAM file (e.g. gs://ref/hg38.fasta)."
    fastafai: "[optional] The index for the FASTA file (default: '<reference_fasta>.fai') (e.g. gs://ref/hg38.fasta.fai)."
    blocklist: "Path to TSV file containing regions where SVs should not be called (see Delly GitHub repo). Format: 'contig | start | end | interval-type' OR only 'contig'."
    sv_type: "SV event type for Delly to search for (e.g. 'DEL', 'DUP', 'INV')."
  }

  Boolean is_bam = basename(crbam, ".bam") + ".bam" == basename(crbam)
  String index_path = if is_bam then crbam + ".bai" else crbam + ".crai"
  File crbam_index_file = select_first([crbami, index_path])
  File ref_index = select_first([fastafai, fasta + ".fai"])

  # ensure there's sufficient disk space
  Float disk_overhead = 10.0
  Float crbam_size = size(crbam, "GiB")
  Float crbam_index_size = size(crbam_index_file, "GiB")
  Float ref_size = size(fasta, "GiB")
  Float ref_index_size = size(ref_index, "GiB")
  Float blocklist_size = size(blocklist, "GiB")
  Int vm_disk_size = ceil(crbam_size + crbam_index_size + ref_size + ref_index_size + blocklist_size + disk_overhead)

  # ensure there's sufficient memory
  Float mem_per_crbam_size = 0.09
  Float mem_crbam_offset = 4.9
  Float mem_size_gb = mem_per_crbam_size * crbam_size + mem_crbam_offset

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_size_gb,
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bcf = "${sample_id}.delly.${sv_type}.bcf"
  }

  command <<<
    set -Eeuo pipefail

    BCF="~{sample_id}.delly.~{sv_type}.bcf"
    echo "Running delly on sv_type=~{sv_type}, output=$BCF"
    delly call \
      -t ~{sv_type} \
      -g "~{fasta}" \
      -x "~{blocklist}" \
      -o "$BCF" \
      -n \
      "~{crbam}"
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

task GatherBCFs {
  input {
    Array[File] bcfs
    String sample_id
    String docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bcfs: "Array of BCF files output by Delly (one per SV event type)."
    sample_id: "Sample name, used for output (e.g. FOO will give a VCF output FOO.vcf.gz)"
  }

  File first_bcf = bcfs[0]
  Int num_bcfs = length(bcfs)

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.7,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "${sample_id}.delly.vcf.gz"
    File index = "${sample_id}.delly.vcf.gz.tbi"
  }

  command <<<
    set -Eeuo pipefail

    BCFS="~{sep='\n' bcfs}"
    echo "Extracting VCF header from ~{first_bcf}"
    bcftools view ~{first_bcf} | grep "^#" > header.vcf
    VCFS="header.vcf"
    for ((BCF_IND=1; BCF_IND <= ~{num_bcfs}; BCF_IND++)); do
      echo "Processing BCF ${BCF_IND}"
      BCF=$(echo "${BCFS}" | sed "$BCF_IND""q;d")
      VCF=$(echo "${BCF}" | sed -e 's/.bcf$/.vcf/')
      echo "Converting $BCF to $VCF, skipping header"
      bcftools view  $BCF | grep -v "^#" > $VCF
      VCFS="$VCFS $VCF"
    done

    VCF_OUT="~{sample_id}.delly.vcf.gz"
    echo "Concatenating VCFs into ${VCF_OUT}"
    cat $VCFS \
      | vcf-sort -c \
      | bcftools reheader -s <(echo "~{sample_id}") \
      | bgzip -c \
      > $VCF_OUT
    echo "Indexing $VCF_OUT"
    tabix "$VCF_OUT"
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
