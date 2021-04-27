#!/usr/bin/env bash

# Build MELT and upload to Google Container Registry

TOOL="melt"
VERSION="2.2.2"
REPO="gcr.io/cpg-common/sv"
TAG="${REPO}/${TOOL}:${VERSION}"

# Build time: 2min
# Needs local MELT tarball to be uploaded.
gcloud builds submit --tag "${TAG}" .
