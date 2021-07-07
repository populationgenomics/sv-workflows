#!/usr/bin/env bash

# Runs PESRCollection workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/analysis-runner-test/PESRCollection",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "78b1d9cc73796631da4522099b5d89bebe5c9ce6",
        "inputs_dict": {
          "PESRCollection.sample_id": "NA12878_nygc",
          "PESRCollection.bam_or_cram_file": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam",
          "PESRCollection.bam_or_cram_index": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam.bai"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": ["../../../gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/Module00a/PESRCollection",
        "description": "PESRCollection on NA12878"
    }'
