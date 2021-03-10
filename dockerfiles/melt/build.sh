#!/usr/bin/env bash

TOOL="melt"
VERSION="2.2.2"
REPO="sv"
PROJECT_ID=$(gcloud config get-value project)
LOCATION="australia-southeast1"
TAG="${LOCATION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/${TOOL}:${VERSION}"

# takes ~2min - needs local MELT tarball to be uploaded
gcloud builds submit --tag "${TAG}" .
