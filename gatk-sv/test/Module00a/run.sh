#!/usr/bin/env bash

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/analysis-runner-test/Module00a",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "2f210a38461ac1cf038f33da881732cc6387d866",
        "inputs_dict": {
            "Module00a.bam_or_cram_file": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam",
            "Module00a.bam_or_cram_index": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam.bai",
            "Module00a.sample_id": "NA12878_nygc"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": [],
        "cwd": "gatk-sv/test/Module00a",
        "description": "Module00a on NA12878"
    }'
