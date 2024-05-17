#!/usr/bin/env bash

# Moves catalog for genotyping into cpg commons bucket so that it can be accessed for genotyping runs from other buckets.

set -ex

gsutil -m cp -r gs://cpg-bioheart-main-analysis/str/polymorphic_run/catalog gs://cpg-common-main/references/str/tob_bioheart_polymorphic_catalog
