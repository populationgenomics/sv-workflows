require(tidyverse)
require(jsonlite, include.only = "read_json")
require(glue, include.only = "glue")

# read json with docker images
gatk_sv <- list(repo = "https://raw.githubusercontent.com/populationgenomics/gatk-sv",
                tag = "v0.16-beta")
gatk_sv_json <-
  glue("{gatk_sv$repo}/{gatk_sv$tag}/input_values/dockers.json") |>
  read_json()

au_ar <- list(loc = "australia-southeast1",
              proj = "fewgenomes",
              repo = "sv")
au_artifact_registry <- glue("{au_ar$loc}-docker.pkg.dev/{au_ar$proj}/{au_ar$repo}")

d <- gatk_sv_json |>
  enframe(name = "image_name", value = "us_gcr") |>
  unnest(us_gcr) |>
  filter(!image_name %in% c("name", "melt_docker")) |>
  mutate(bname = basename(us_gcr),
         au_ar = glue("{au_artifact_registry}/{bname}")) |>
  select(-bname)

# for (i in seq_len(nrow(d))) {
for (i in c(9, 10, 11, 12, 13, 14, 15)) {
  cat(glue("pulling {d$us_gcr[i]}"), "\n")
  system(glue("docker pull {d$us_gcr[i]}"))
  system(glue("docker tag {d$us_gcr[i]} {d$au_ar[i]}"))
  cat(glue("pushing {d$au_ar[i]}"), "\n")
  system(glue("docker push {d$au_ar[i]}"))
}
