#!/usr/bin/env bash

# Moves Hail Matrix Table containing STR genotypes into test bucket

set -ex

gsutil -m cp -r gs://cpg-bioheart-main/str/polymorphic_run_n2045/mt/v1/str.mt gs://cpg-bioheart-test/str/polymorphic_run_n2045/mt/v1/str.mt
