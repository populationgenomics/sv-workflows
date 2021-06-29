#!/usr/bin/env bash

set -euo pipefail

# Build MELT and upload to Google Container Registry

REGISTRY="gcr.io/cpg-common/sv"
TOOL="melt"
VERSION="2.2.2"
TAG="${REGISTRY}/${TOOL}:${VERSION}"
TARBALL="MELTv${VERSION}.tar.gz"

# Needs local MELT tarball
if [[ ! -f  "$TARBALL" ]]; then
    echo "Downloading MELT tarball from Google Cloud Storage"
    gsutil cp gs://cpg-reference/sv/MELT/${TARBALL} .
fi

docker image build --tag "${TAG}" --build-arg MELT_RELEASE="${VERSION}" .
docker image push "${TAG}"
