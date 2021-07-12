version 1.0

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/MakeBincovMatrix.wdl" as mbm

workflow MakeBincovMatrix {
  input {
    Array[String] samples
    Array[File] count_files
    Array[String]? bincov_matrix_samples
    File? bincov_matrix
    String batch
    Int? binsize
    Int? disk_overhead_gb
    String sv_base_mini_docker
    String sv_base_docker
    RuntimeAttr? runtime_attr_override
  }

  call mbm.SetBins {
    input:
      count_file = all_count_files[0],
      binsize = binsize,
      bincov_matrix_samples = bincov_matrix_samples,
      sv_base_mini_docker = sv_base_mini_docker,
      disk_overhead_gb = disk_overhead_gb,
      runtime_attr_override = runtime_attr_override
  }
  if(defined(bincov_matrix_samples)) {
    String bincov_matrix_header = read_lines(SetBins.bincov_matrix_header_file)[0]
  }

  Array[String]+ all_samples = flatten([samples, select_all([bincov_matrix_header])])
  Array[File]+ all_count_files = flatten([count_files, select_all([bincov_matrix])])

  scatter(i in range(length(all_count_files))) {
    call mbm.MakeBincovMatrixColumns {
      input:
        count_file = all_count_files[i],
        sample = all_samples[i],
        binsize = SetBins.out_binsize,
        bin_locs = SetBins.bin_locs,
        disk_overhead_gb = disk_overhead_gb,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  call mbm.ZPaste {
    input:
      column_files = flatten([[SetBins.bin_locs], MakeBincovMatrixColumns.bincov_bed]),
      matrix_file_name = "~{batch}.RD.txt.gz",
      disk_overhead_gb = disk_overhead_gb,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File merged_bincov = ZPaste.matrix_file
    File merged_bincov_idx = ZPaste.matrix_file_idx
  }
}
