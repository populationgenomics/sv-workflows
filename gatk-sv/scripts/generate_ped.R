require(tidyverse)
require(jsonlite)
require(here)

cpg2tob <- here("nogit/cpg2tob_96samples_idmap.json") |>
  read_json() |>
  enframe(name = "cpgid", value = "tobid") |>
  unnest(tobid)

reported_sex <- here("nogit/reported_sex.tsv") |>
  read_tsv(skip = 2,
           col_types = cols(.default = "c"),
           col_names = c("sample_id", "sample.sample_name", "gender",
                         paste0("nanodrop_", c("conc", "260/280", "230/280")),
                         paste0("qubit_", c("assay_conc", "stock_conc")),
                         "volume", "location")) %>%
  filter(!is.na(sample_id)) %>%
  select(tobid = sample.sample_name, reported_sex = gender)

# from GATK-SV's EvidenceQC workflow
gatksv_sex <- here("nogit/sample_sex_assignments.txt.gz") |>
  read_tsv(col_types = cols(.default = "c")) |>
  select(cpgid = `#ID`, sex = Assignment)

d <- cpg2tob |>
  left_join(reported_sex, by = "tobid") |>
  left_join(gatksv_sex, by = "cpgid")

# known sex labeling error for TOB1856
d |>
  filter(reported_sex == "F", sex == "MALE")

# Create pedigree
d |>
  mutate(FamID = cpgid,
         IndID = cpgid,
         PatID = 0,
         MatID = 0,
         Sex = if_else(sex == "MALE", 1, if_else(sex == "FEMALE", 2, 0)),
         Pheno = 0) |>
  select(FamID:Pheno) |>
  write_tsv(file = here("nogit/tob96.ped"), col_names = FALSE)
