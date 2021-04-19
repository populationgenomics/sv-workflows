#!/usr/bin/env bash

# Pulling and pushing based on tagged versions in
# https://github.com/broadinstitute/gatk-sv/blob/master/input_values/dockers.json

TOOL="sv-base-mini"
VERSION="rlc_posthoc_filtering_cnv_mcnv_compatability_9a8561"
AU_REPO="sv"
US_REPO="gatk-sv"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="australia-southeast1"

US_TAG="us.gcr.io/broad-dsde-methods/${US_REPO}/${TOOL}:${VERSION}"
AU_TAG="${LOCATION}-docker.pkg.dev/${PROJECT_ID}/${AU_REPO}/${TOOL}:${VERSION}"


docker pull $US_TAG
docker tag $US_TAG $AU_TAG
docker push $AU_TAG
