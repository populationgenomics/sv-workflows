#!/usr/bin/env Rscript

# skopeo required for copying images
system("micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge skopeo")

require(dplyr)
require(tidyr)
require(tibble)
require(jsonlite, include.only = "read_json")
require(glue, include.only = "glue")

# read json with docker images
gatk_sv <- list(repo = "https://raw.githubusercontent.com/populationgenomics/gatk-sv",
                tag = "v0.18-beta")
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
  filter(!image_name %in% c("name", "melt_docker", "cloud_sdk_docker", "linux_docker")) %>%
  mutate(bname = basename(us_gcr)) %>%
  select(us_gcr, bname) %>%
  distinct()

# copy to AU AR
system("gcloud auth configure-docker australia-southeast1-docker.pkg.dev")
for (i in seq_len(nrow(d))) {
  system(glue("skopeo copy docker://{d$us_gcr[i]} docker://{au_artifact_registry}/{d$bname[i]}"))
}
