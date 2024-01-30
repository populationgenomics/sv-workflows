#!/usr/bin/env bash

# Moves catalog for genotyping into TOB bucket so that it can be accessed for the TOB run

set -ex

gsutil -m cp -r gs://cpg-bioheart-main-analysis/str/polymorphic_run/catalog gs://cpg-tob-wgs-main-analysis/str/polymorphic_run/catalog
