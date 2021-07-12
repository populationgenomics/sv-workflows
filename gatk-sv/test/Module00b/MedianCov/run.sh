#!/usr/bin/env bash

# Runs MakeBincovMatrix workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/MedianCov",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "8f9b97a4c8a88b75d655cefb421e4ffb2ee9ff2b",
        "inputs_dict": {
          "MedianCov.cohort_id": "test_NA12878"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00b/MedianCov",
        "description": "MedianCov on NA12878"
    }'
