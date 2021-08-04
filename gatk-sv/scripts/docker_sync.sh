#!/usr/bin/env bash

micromamba install -y --prefix $MAMBA_ROOT_PREFIX -y -c conda-forge \
    skopeo

skopeo inspect docker://registry.fedoraproject.org/fedora:latest
