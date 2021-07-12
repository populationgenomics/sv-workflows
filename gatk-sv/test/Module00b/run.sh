#!/usr/bin/env bash

# Runs Module00b workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/Module00b",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "fcb3ad57515df466552a28c15acd9b6755406cd0",
        "inputs_dict": {
          "Module00b.samples": ["NA12878"],
          "Module00b.batch": "NA12878",
          "Module00b.counts": ["gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/gatk-sv/CollectCounts/NA12878.counts.tsv.gz"]
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00b",
        "description": "Module00b on NA12878"
    }'

