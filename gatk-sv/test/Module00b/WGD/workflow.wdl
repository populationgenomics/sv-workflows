version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.15.1-beta/wdl/Structs.wdl"
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.15.1-beta/wdl/WGD.wdl" as wgd

workflow WGD {
  input {
    String batch
    File wgd_scoring_mask
    File bincov_matrix
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_build
    RuntimeAttr? runtime_attr_score
  }

  call wgd.BuildWGDMatrix {
    input:
      bincov_matrix = bincov_matrix,
      wgd_scoring_mask = wgd_scoring_mask,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      batch = batch
  }

  call wgd.WGDScore {
    input:
      wgd_scoring_mask = wgd_scoring_mask,
      WGD_matrix = BuildWGDMatrix.WGD_matrix,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      batch = batch
  }

  output {
    File WGD_dist = WGDScore.WGD_dist
    File WGD_matrix = BuildWGDMatrix.WGD_matrix
    File WGD_scores = WGDScore.WGD_scores
  }
}
