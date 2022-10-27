#!/usr/bin/env bash

# Runs GATK-SV SingleSample Pipeline using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/gatksv_single/NA12878",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "4e554f59eafc7f7b220decdce177649d9bf400d4",
        "inputs_dict": {
          "GATKSVPipelineSingleSample.sample_id" : "NA12878",
          "GATKSVPipelineSingleSample.batch" : "test_NA12878",
          "GATKSVPipelineSingleSample.bam_or_cram_file" : "gs://cpg-fewgenomes-test/pdiakumis/data/NA12878/NA12878.final.bam"
        },
        "input_json_paths": ["inputs/general.json"],
        "workflow": "GATKSVPipelineSingleSample.wdl",
        "dependencies": ["gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test/GATKSVPipelineSingleSample",
        "description": "Single sample pipeline on NA12878"
    }'
