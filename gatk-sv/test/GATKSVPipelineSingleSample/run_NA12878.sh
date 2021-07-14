#!/usr/bin/env bash

# Runs GATK-SV SingleSample Pipeline using analysis-runner

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/test1/gatksv_single",
        "dataset": "fewgenomes",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "9fb01edc3bf6dfb234fc522f61ccf9167d3a2f9a",
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

