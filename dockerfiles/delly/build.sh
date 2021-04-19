#!/usr/bin/env bash

# Build and upload to Google Cloud Container Registry

TOOL="delly"
VERSION="0.8.7"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="us"
TAG="${LOCATION}.gcr.io/${PROJECT_ID}/${TOOL}:${VERSION}"

# Build time: 4min
gcloud builds submit --tag "${TAG}"
