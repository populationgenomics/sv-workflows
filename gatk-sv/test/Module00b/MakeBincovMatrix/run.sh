#!/usr/bin/env bash

# Runs MakeBincovMatrix workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/MakeBincovMatrix",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "fcece36a1cd811efbcb03b6220d465580111a4c6",
        "inputs_dict": {
          "MakeBincovMatrix.samples": ["NA12878"],
          "MakeBincovMatrix.count_files": ["gs://cpg-fewgenomes-test/pdiakumis/test1/Module00a/Module00a/74ed65a0-245b-4db5-a2b0-b8a23c719bc1/call-CollectCounts/NA12878_test1.counts.tsv.gz"],
          "MakeBincovMatrix.batch": "test_NA12878"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00b/MakeBincovMatrix",
        "description": "MakeBincovMatrix on NA12878"
    }'
