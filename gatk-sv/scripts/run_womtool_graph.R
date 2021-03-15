require(tidyverse)
require(here)
require(glue)

# grab initial WDL files to plot WOMTool graphs
wdl_files <-
  tibble::tibble(
    basename = c(
      "Module00a.wdl", "Module00b.wdl", "Module00c.wdl", "Module01.wdl",
      "CramToBam.wdl", "PESRCollection.wdl",
      "Whamg.wdl", "Delly.wdl", "Manta.wdl", "MELT.wdl")) %>%
  dplyr::mutate(
    fullname = file.path(here("gatk-sv/scripts/wdl"), basename),
    figure = file.path(here("gatk-sv/figures"), paste0(basename, ".graph.svg")))

conda_env <- "/Users/peterd/conda/envs/wdl/bin"
for (i in seq_len(nrow(wdl_files))) {
  system(glue::glue("{conda_env}/womtool graph {wdl_files$fullname[i]} | ",
                    "dot -Tsvg -o {wdl_files$figure[i]}"))
}
