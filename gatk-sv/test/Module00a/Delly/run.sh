#!/usr/bin/env bash

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/analysis-runner-test/delly",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "423da419dd78571eda91f7c4fa2729bcd20b49e9",
        "inputs_dict": {
          "Delly.sample_id": "NA12878_nygc",
          "Delly.bam_or_cram_file": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam",
          "Delly.bam_or_cram_index": "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam.bai"
        },
        "input_json_paths": ["inputs_general.json"],
        "workflow": "workflow.wdl",
        "dependencies": ["../../../gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/Module00a/Delly",
        "description": "Delly on NA12878"
    }'

