require(tidyverse)
require(here)
require(glue)

conda_env <- "wdl"
womtool <- glue("{Sys.getenv('CONDA_PREFIX')}/envs/{conda_env}/bin/womtool")

# grab initial WDL files to plot WOMTool graphs
wdl_files <-
  tibble(
    basename = c(
      "Module00a.wdl", "Module00b.wdl", "Module00c.wdl", "Module01.wdl",
      "CramToBam.wdl", "PESRCollection.wdl",
      "Whamg.wdl", "Manta.wdl", "CNMOPS.wdl")) %>%
  mutate(
    fullname = glue("{here('gatk-sv/gatk-sv-git/wdl')}/{basename}"),
    figure = glue("{here('gatk-sv/figures')}/{basename}.graph.svg"))

for (i in seq_len(nrow(wdl_files))) {
  system(glue("{womtool} graph {wdl_files$fullname[i]} | ",
              "dot -Tsvg -o {wdl_files$figure[i]}"))
}
