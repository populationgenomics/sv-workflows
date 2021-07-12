#!/usr/bin/env bash

# Runs WGD workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/WGD",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "d59506f47e179ac4f4b3707f4d065f0e439e8281",
        "inputs_dict": {
          "WGD.batch": "test_NA12878",
          "WGD.bincov_matrix": "gs://cpg-fewgenomes-test/pdiakumis/test1/MakeBincovMatrix/MakeBincovMatrix/3f6f638a-93c5-4bc2-81e2-1429f5c6fa28/call-ZPaste/test_NA12878.RD.txt.gz"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00b/WGD",
        "description": "WGD on NA12878"
    }'

