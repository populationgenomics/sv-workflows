##  This WDL implements workflow for ExpansionHunter.

version 1.0

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.17.1-beta/wdl/Structs.wdl"

struct FilenamePostfixes {
    String locus
    String motif
    String profile
    String merged_profile
    Int profile_len
}

workflow ExpansionHunter {

    input {
        File bam_or_cram
        File? bam_or_cram_index
        File reference_fasta
        File? reference_fasta_index
        File variant_catalog
        String sex
        String output_prefix
        String docker_file
        RuntimeAttr? runtime_attr
    }

    Boolean is_bam = basename(bam_or_cram, ".bam") + ".bam" == basename(bam_or_cram)
    File bam_or_cram_index_ =
        if defined(bam_or_cram_index) then
            select_first([bam_or_cram_index])
        else
            bam_or_cram + if is_bam then ".bai" else ".crai"

    File reference_fasta_index_ = select_first([
        reference_fasta_index,
        reference_fasta + ".fai"])

    call RunExpansionHunter {
        input:
            bam_or_cram = bam_or_cram,
            bam_or_cram_index = bam_or_cram_index_,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index_,
            variant_catalog = variant_catalog,
            sex = sex,
            output_prefix = output_prefix,
            docker_file = docker_file,
            runtime_attr_override = runtime_attr,
    }

    output {
        File json = RunExpansionHunter.json
        File vcf = RunExpansionHunter.vcf
        File overlapping_reads = RunExpansionHunter.overlapping_reads
    }
}

task RunExpansionHunter {
    input {
        File bam_or_cram
        File bam_or_cram_index
        File reference_fasta
        File reference_fasta_index
        File variant_catalog
        String sex
        String docker_file
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }


    output {
        File json = "${output_prefix}.json"
        File vcf = "${output_prefix}.vcf"
        File overlapping_reads = "${output_prefix}_realigned.bam"
    }

    command <<<
        set -euxo pipefail

        ExpansionHunter \
            --reads ~{bam_or_cram} \
            --reference ~{reference_fasta} \
            --variant-catalog ~{variant_catalog} \
            --output-prefix ~{output_prefix} \
            --sex ~{sex}
    >>>

    RuntimeAttr runtime_attr_str_profile_default = object {
        cpu_cores: 1,
        mem_gb: 32,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1,
        disk_gb: 10 + ceil(size([
            bam_or_cram,
            bam_or_cram_index,
            reference_fasta,
            reference_fasta_index], "GiB"))
    }
    RuntimeAttr runtime_attr = select_first([
        runtime_attr_override,
        runtime_attr_str_profile_default])

    runtime {
        docker: docker_file
        cpu: runtime_attr.cpu_cores
        memory: runtime_attr.mem_gb + " GiB"
        disks: "local-disk " + runtime_attr.disk_gb + " HDD"
        bootDiskSizeGb: runtime_attr.boot_disk_gb
        preemptible: runtime_attr.preemptible_tries
        maxRetries: runtime_attr.max_retries
    }
}
