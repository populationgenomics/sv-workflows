#!/usr/bin/env bash

set -euo pipefail

# Build MELT and upload to Google Artifact Registry

LOCATION="australia-southeast1"
PROJECT="peter-dev-302805"
TOOL="melt"
VERSION="2.2.2"
TARBALL="MELTv${VERSION}.tar.gz"
REPO="test"
TAG="${LOCATION}-docker.pkg.dev/${PROJECT}/${REPO}/${TOOL}:${VERSION}"

# Needs local MELT tarball
if [[ ! -f  "$TARBALL" ]]; then
    echo "Downloading MELT tarball from Google Cloud Storage"
    gsutil cp gs://cpg-common-main/references/sv/MELT/${TARBALL} .
fi

docker image build --tag "${TAG}" --build-arg MELT_RELEASE="${VERSION}" .
docker image push "${TAG}"
