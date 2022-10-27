#!/usr/bin/env bash

# Runs Manta workflow using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/manta",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "141ce634a3e55f6f8b6b1f24a72b208c35e9ab7d",
        "inputs_dict": {
          "Manta.sample_id": "NA12878_test1",
          "Manta.bam_or_cram_file": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam",
          "Manta.bam_or_cram_index": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam.bai"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00a/Manta",
        "description": "Manta on NA12878"
    }'
