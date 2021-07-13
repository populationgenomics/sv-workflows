#!/usr/bin/env bash

# Runs Module00b workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/gatksv_single",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "d27963f0a31d9789173793eb41eac780ec49810b",
        "input_json_paths": ["inputs/GATKSVPipelineSingleSample.ref_panel_1kg.na12878.no_melt.json"],
        "workflow": "GATKSVPipelineSingleSample.wdl",
        "dependencies": ["gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/GATKSVPipelineSingleSample",
        "description": "Single sample pipeline on NA12878"
    }'

