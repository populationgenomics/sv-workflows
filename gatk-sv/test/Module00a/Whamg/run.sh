#!/usr/bin/env bash

# Runs Whamg workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/analysis-runner-test/whamg",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "554d1481625994c98d3bcddaf17ffbd30baa6e89",
        "inputs_dict": {
          "Whamg.sample_id": "NA12878_nygc",
          "Whamg.bam_or_cram_file": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam",
          "Whamg.bam_or_cram_index": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam.bai"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": ["../../../gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/Module00a/Whamg",
        "description": "Whamg on NA12878"
    }'

