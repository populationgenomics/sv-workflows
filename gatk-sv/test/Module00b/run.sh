#!/usr/bin/env bash

# Runs WGD workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/Module00b",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "3020b9088a187ae036b035619f37e94962b5f3a0",
        "inputs_dict": {
          "Module00b.samples": ["NA12878"],
          "Module00b.batch": "test_NA12878",
          "Module00b.counts": ["gs://cpg-fewgenomes-test/pdiakumis/test1/Module00a/Module00a/74ed65a0-245b-4db5-a2b0-b8a23c719bc1/call-CollectCounts/NA12878_test1.counts.tsv.gz"]
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00b",
        "description": "Module00b on NA12878"
    }'

