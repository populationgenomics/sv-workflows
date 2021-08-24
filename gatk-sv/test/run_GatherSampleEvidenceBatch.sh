#!/usr/bin/env bash

curl --location \
    --request POST 'https://server-a2pko7ameq-ts.a.run.app/cromwell' \
    --header "Authorization: Bearer $(gcloud auth print-identity-token)" \
    --header 'Content-Type: application/json' \
    --data-raw '{
        "output": "pdiakumis/gatksvbatch",
        "dataset": "tob-wgs",
        "repo": "sv-workflows",
        "accessLevel": "test",
        "commit": "256fe8f175b588bb4876b551d5d33a4d814bfcd9",
        "inputs_dict": {
          "GatherSampleEvidenceBatch.sample_ids": [
              "TOB1520",
              "TOB1521",
              "TOB1522"
          ],
          "GatherSampleEvidenceBatch.bam_or_cram_files": [
              "gs://cpg-tob-wgs-test/cram/batch1/TOB1520.cram",
              "gs://cpg-tob-wgs-test/cram/batch1/TOB1521.cram",
              "gs://cpg-tob-wgs-test/cram/batch1/TOB1522.cram"
          ],
          "GatherSampleEvidenceBatch.batch": "TOB_test1"
        },
        "input_json_paths": ["GatherSampleEvidenceBatch_inputs.json"],
        "workflow": "GatherSampleEvidenceBatch.wdl",
        "dependencies": ["gatk-sv-git/wdl"],
        "cwd": "gatk-sv/test",
        "description": "GatherSampleEvidenceBatch"
    }'
