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
              proj = "fewgenomes",
              repo = "sv")
au_artifact_registry <- glue("{au_ar$loc}-docker.pkg.dev/{au_ar$proj}/{au_ar$repo}")

# dataframe linking US GCR image to AU AR copy
d <- gatk_sv_json %>%
  enframe(name = "image_name", value = "us_gcr") %>%
  unnest(us_gcr) %>%
  filter(!image_name %in% c("name", "melt_docker")) %>%
  mutate(bname = basename(us_gcr),
         au_ar = glue("{au_artifact_registry}/{bname}")) %>%
  select(-bname)

# pull from US and push to AU
for (i in seq_len(d)) {
  cat(glue("pulling {d$us_gcr[i]}"), "\n")
  system(glue("docker pull {d$us_gcr[i]}"))
  system(glue("docker tag {d$us_gcr[i]} {d$au_ar[i]}"))
  cat(glue("pushing {d$au_ar[i]}"), "\n")
  system(glue("docker push {d$au_ar[i]}"))
}
