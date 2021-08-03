#!/usr/bin/env Rscript

require(tidyverse)
require(jsonlite, include.only = "read_json")
require(glue, include.only = "glue")

# read json with docker images
gatk_sv <- list(repo = "https://raw.githubusercontent.com/populationgenomics/gatk-sv",
                tag = "v0.16-beta")
gatk_sv_json <-
  glue("{gatk_sv$repo}/{gatk_sv$tag}/input_values/dockers.json") %>%
  read_json()

# AR repo to hold the images
au_ar <- list(loc = "australia-southeast1",
              proj = "cpg-common",
              repo = "images/sv")
au_artifact_registry <- glue("{au_ar$loc}-docker.pkg.dev/{au_ar$proj}/{au_ar$repo}")

# dataframe with US GCR images
d <- gatk_sv_json %>%
  enframe(name = "image_name", value = "us_gcr") %>%
  unnest(us_gcr) %>%
  filter(!image_name %in% c("name", "melt_docker")) %>%
  mutate(bname = basename(us_gcr)) %>%
  select(us_gcr, bname)

# copy to AU AR
DOCKERFILE <- "Dockerfile"
for (i in seq_len(nrow(d))) {
  if (file.exists(DOCKERFILE)) {
    file.remove(DOCKERFILE)
  }
  write_lines(glue("FROM {d$us_gcr[i]}"), file = DOCKERFILE)
  system(glue("gcloud builds submit --tag {au_artifact_registry}/{d$bname[i]} ."))
  file.remove(DOCKERFILE)
}
