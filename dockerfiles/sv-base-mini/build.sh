#!/usr/bin/env bash

# Pulling and pushing based on tagged versions (as of 2021-Apr-19) in
# https://github.com/broadinstitute/gatk-sv/blob/master/input_values/dockers.json

TOOL="sv-base-mini"
VERSION="rlc_posthoc_filtering_cnv_mcnv_compatability_9a8561"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="us"

ORIGINAL_TAG="us.gcr.io/broad-dsde-methods/gatk-sv/${TOOL}:${VERSION}"
NEW_TAG="${LOCATION}.gcr.io/${PROJECT_ID}/${TOOL}:${VERSION}"


docker pull $ORIGINAL_TAG
docker tag $ORIGINAL_TAG $NEW_TAG
docker push $NEW_TAG
