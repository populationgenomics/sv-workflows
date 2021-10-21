#!/usr/bin/env bash

# Runs GATK-SV SingleSample Pipeline using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "vladsavelyev/gatksvbatch_clusterbatch",
        "dataset": "tob-wgs",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "52bb628a7739b21aed7c3ec8c0d9a2b586d5f3b6",
        "input_json_paths": ["inputs/ClusterBatch.json"],
        "workflow": "gatk-sv-git/wdl/ClusterBatch.wdl",
        "dependencies": ["gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/ClusterBatch",
        "description": "ClusterBatch for TOB-WGS"
    }'
