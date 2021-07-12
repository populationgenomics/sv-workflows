version 1.0

# Workflow definition for Calculating Median Coverage 

import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/populationgenomics/gatk-sv/v0.15.1-beta/wdl/MedianCov.wdl" as mc

workflow MedianCov {
  input {
    File bincov_matrix
    String cohort_id
    Float? mem_gb_override
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr
  }

  call mc.CalcMedCov {
    input:
      bincov_matrix = bincov_matrix,
      cohort_id = cohort_id,
      mem_gb_override = mem_gb_override,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_override = runtime_attr
  }

  output {
    File medianCov = CalcMedCov.medianCov
  }
}
