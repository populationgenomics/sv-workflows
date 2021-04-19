#!/usr/bin/env bash

# Build and upload to Google Cloud Container Registry

TOOL="wham"
VERSION="1.8.0"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="us"
TAG="${LOCATION}.gcr.io/${PROJECT_ID}/${TOOL}:${VERSION}"

# Build time: 1min
# Needs local whamg binary to be uploaded (grabbed from GATK-SV repo).
gcloud builds submit --tag "${TAG}" .
