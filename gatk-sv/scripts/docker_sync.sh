#!/usr/bin/env bash

set -ex

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge skopeo

skopeo copy docker://docker.io/library/alpine:3.12 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/sv/alpine:foo
