#!/usr/bin/env bash

# Copies Hail Matrix Table containing STR genotypes into test bucket

set -ex

#gsutil -m cp -r gs://cpg-bioheart-main/str/polymorphic_run_n2045/mt/v1/str.mt gs://cpg-bioheart-test/str/polymorphic_run_n2045/mt/v1/str.mt

# copy cohort specific mts (changing approach to analyse each cohort separately because of distinct batch effects)

gcloud storage cp -r gs://cpg-bioheart-main/str/polymorphic_run_n990_bioheart_only/annotated_mt/v1/str_annotated.mt gs://cpg-bioheart-test/str/polymorphic_run_n990_bioheart_only/annotated_mt/v1/str_annotated.mt
gcloud storage cp -r gs://cpg-bioheart-main/str/polymorphic_run_n1055_tob_only/annotated_mt/v1/str_annotated.mt gs://cpg-bioheart-test/str/polymorphic_run_n1055_tob_only/annotated_mt/v1/str_annotated.mt
