#!/usr/bin/env bash

set -ex

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://weisburd/gangstr:v2.5 docker://australia-southeast1-docker.pkg.dev/cpg-common/images/gangstr:v2.5