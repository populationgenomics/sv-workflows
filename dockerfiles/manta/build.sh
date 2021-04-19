#!/usr/bin/env bash

TOOL="manta"
VERSION="1.6.0"
REPO="sv"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="australia-southeast1"
TAG="${LOCATION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/${TOOL}:${VERSION}"

# takes ~3min
gcloud builds submit --tag "${TAG}"
