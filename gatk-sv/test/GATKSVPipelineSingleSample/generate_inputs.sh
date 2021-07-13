#!/usr/bin/env bash

GATK_SV_REPO="$(PWD)/gatk-sv-git"
python ${GATK_SV_REPO}/scripts/inputs/build_inputs.py \
       ${GATK_SV_REPO}/input_values \
       ${GATK_SV_REPO}/input_templates \
       $(PWD)/inputs2 \
       -a '{ "single_sample" : "test_single_sample_NA12878", "ref_panel" : "ref_panel_1kg_v2"}'
