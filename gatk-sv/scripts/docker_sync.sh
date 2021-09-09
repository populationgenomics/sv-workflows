#!/usr/bin/env bash

set -ex

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge skopeo

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://australia-southeast1-docker.pkg.dev/peter-dev-302805/test/melt:2.2.2 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/sv/melt:2.2.2
