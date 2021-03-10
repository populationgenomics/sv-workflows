#!/usr/bin/env bash

TOOL="delly"
VERSION="0.8.7"
REPO="sv"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="australia-southeast1"
TAG="${LOCATION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/${TOOL}:${VERSION}"

# takes ~4min
gcloud builds submit --tag "${TAG}"
