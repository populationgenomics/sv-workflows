#!/usr/bin/env bash

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge \
    skopeo

echo ${GOOGLE_APPLICATION_CREDENTIALS}
#skopeo inspect docker://australia-southeast1-docker.pkg.dev/peter-dev-302805/tmp/ubuntu1804:foo
skopeo copy docker://docker.io/library/alpine:3.12 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/sv/alpine:3.12
