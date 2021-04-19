#!/usr/bin/env bash

# Build and upload to Google Cloud Container Registry

TOOL="manta"
VERSION="1.6.0"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="us"
TAG="${LOCATION}.gcr.io/${PROJECT_ID}/${TOOL}:${VERSION}"

# Build time: 3min
gcloud builds submit --tag "${TAG}"
