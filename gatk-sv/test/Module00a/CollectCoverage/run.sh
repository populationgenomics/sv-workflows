#!/usr/bin/env bash

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/analysis-runner-test/collectCoverage",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "caffefcf779ce85d4c6d1d20655a59d69510a9d6",
        "inputs_dict": {
          "CollectCoverage.sample_id": "NA12878_nygc",
          "CollectCoverage.bam_or_cram_file": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam",
          "CollectCoverage.bam_or_cram_index": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam.bai"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": ["../../../gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/Module00a/CollectCoverage",
        "description": "CollectCoverage on NA12878"
    }'
