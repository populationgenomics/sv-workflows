#!/usr/bin/env bash

# Build and upload to Google Cloud Container Registry

TOOL="melt"
VERSION="2.2.2"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="us"
TAG="${LOCATION}.gcr.io/${PROJECT_ID}/${TOOL}:${VERSION}"

# Build time: 2min
# Needs local MELT tarball to be uploaded.
gcloud builds submit --tag "${TAG}" .
