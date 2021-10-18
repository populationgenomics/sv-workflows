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
        "commit": "4e554f59eafc7f7b220decdce177649d9bf400d4",
        "input_json_paths": ["inputs/ClusterBatch.json"],
        "workflow": "gatk-sv-git/wdl/ClusterBatch.wdl",
        "dependencies": ["gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/ClusterBatch",
        "description": "ClusterBatch for TOB-WGS"
    }'
