#!/usr/bin/env bash

# Runs MedianCov workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/MedianCov",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "fcb3ad57515df466552a28c15acd9b6755406cd0",
        "inputs_dict": {
          "MedianCov.cohort_id": "test_NA12878"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00b/MedianCov",
        "description": "MedianCov on NA12878"
    }'
