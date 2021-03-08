# Dockerfiles

We use Google
[Cloud Build](https://cloud.google.com/build/docs/quickstart-build) to build all
Docker images and upload to Google
[Artifact Registry](https://cloud.google.com/artifact-registry/docs/docker/quickstart).

We build using the following:

```shell
gcloud builds submit --config cloudbuild.yaml
```

A typical `cloudbuild.yaml` is (assuming a `gatk-sv` Artifact Registry
repository has been created under the default Google project beforehand):

```yaml
steps:
  - name: "gcr.io/cloud-builders/docker"
    args:
      [
        "build",
        "-t",
        "australia-southeast1-docker.pkg.dev/$PROJECT_ID/gatk-sv/IMG:TAG",
        ".",
      ]
images:
  - "australia-southeast1-docker.pkg.dev/$PROJECT_ID/gatk-sv/IMG:TAG"
```
