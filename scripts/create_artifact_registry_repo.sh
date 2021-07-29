#!/usr/bin/env bash

LOCATION="australia-southeast1"
PROJECT="fewgenomes"
REPO="sv"

gcloud artifacts repositories create ${REPO} \
    --repository-format "docker" \
    --location ${LOCATION} \
    --description "Structural Variation" \
    --project ${PROJECT}
