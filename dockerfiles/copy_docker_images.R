#---- Create Dockerfile JSON ----#
require(tidyverse)
require(glue)
require(here)

#---- Preparation ----#
cpg_registry <- "gcr.io/peter-dev-302805/sv"
dockers <- tibble::tribble(
  ~name,          ~version,                                              ~registry,
  "delly",        "0.8.7--he03298f_1",                                   "quay.io/biocontainers",
  "manta",        "1.6.0--h9ee0642_1",                                   "quay.io/biocontainers",
  "wham",         "1.8.0.1.2017.05.03--h8b12597_1",                      "quay.io/biocontainers",
  "sv-base-mini", "rlc_posthoc_filtering_cnv_mcnv_compatability_9a8561", "us.gcr.io/broad-dsde-methods/gatk-sv",
) %>%
  dplyr::mutate(
    source = glue("{registry}/{name}:{version}"),
    target = glue("{cpg_registry}/{name}:{version}")
  ) %>%
  dplyr::select(name, source, target)

# perform following linearly (doesn't take that long)
walk(glue_data(dockers, "docker image pull {source}"), system)
walk(glue_data(dockers, "docker image tag {source} {target}"), system)
walk(glue_data(dockers, "docker image push {target}"), system)

readr::write_tsv(dockers, here::here("dockerfiles/dockers.tsv"))
