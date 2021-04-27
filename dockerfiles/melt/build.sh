#!/usr/bin/env bash

set -euo pipefail

# Build MELT and upload to Google Container Registry

TOOL="melt"
VERSION="2.2.2"
REGISTRY="gcr.io/cpg-common/sv"
TAG="${REGISTRY}/${TOOL}:${VERSION}"
TAR="MELTv${VERSION}.tar.gz"

if [[ ! -f  "$TAR" ]]; then
    echo "Downloading MELT tarball from Google Cloud Storage"
    gsutil cp gs://cpg-sv-workflows/MELT/${TAR} .
fi

# Build time: 2min
# Needs local MELT tarball to be uploaded.
gcloud builds submit --tag "${TAG}" .
