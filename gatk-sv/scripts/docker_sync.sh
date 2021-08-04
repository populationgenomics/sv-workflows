#!/usr/bin/env bash

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge \
    skopeo

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://docker.io/library/alpine:3.12 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/sv/alpine:3.12
