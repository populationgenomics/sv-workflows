#!/usr/bin/env bash

PROJECT_ID=$( gcloud config get-value project )
LOCATION="australia-southeast1"
REPO="sv"

IMG="delly"
TAG="0.8.7"

# takes ~4min
gcloud builds submit \
    --tag ${LOCATION}-docker.pkg.dev/$PROJECT_ID/${REPO}/${IMG}:${TAG}
