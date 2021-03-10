# Dockerfiles

We use Google
[Cloud Build](https://cloud.google.com/build/docs/quickstart-build) to build all
Docker images and upload to Google
[Artifact Registry](https://cloud.google.com/artifact-registry/docs/docker/quickstart).

We build using the following (assuming the Dockerfile is in the current working directory):

- create AR repo (done once)

```shell
gcloud artifacts repositories create "sv" \
    --repository-format "docker" \
    --location "australia-southeast1" \
    --description "SV Dockers" \
    --project "<PROJECT_ID>"
```

- build Docker image and upload to specified AR repo under given tag

```shell
gcloud builds submit --tag australia-southeast1-docker.pkg.dev/PROJECT_ID/sv/IMG:TAG
```
