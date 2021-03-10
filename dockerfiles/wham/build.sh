#!/usr/bin/env bash

TOOL="wham"
VERSION="1.8.0"
REPO="sv"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="australia-southeast1"
TAG="${LOCATION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/${TOOL}:${VERSION}"

# takes ~1min - needs local whamg binary to be uploaded (grabbed from GATK-SV)
gcloud builds submit --tag "${TAG}" .

