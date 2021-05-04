version 1.0

import "../delly.wdl" as delly

workflow RunGatherBCFs {

  input {
    Array[File] bcfs
    String sample_id
    String docker
  }


  call delly.GatherBCFs {
      input:
        bcfs = bcfs,
        sample_id = sample_id,
        docker = docker
      }
}
