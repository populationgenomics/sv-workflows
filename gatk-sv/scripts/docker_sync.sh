#!/usr/bin/env bash

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge \
    skopeo

skopeo inspect docker://australia-southeast1-docker.pkg.dev/peter-dev-302805/tmp/ubuntu1804:foo
