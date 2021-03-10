# Dockerfiles

We use [Cloud Build](https://cloud.google.com/build/docs/quickstart-build) to
build all Docker images for the main structural variant tools; these then get
automatically uploaded to
[Artifact Registry](https://cloud.google.com/artifact-registry/docs/docker/quickstart)
(AR).

- create AR repo named `sv` under a given Google Cloud project (done once):

```shell
gcloud artifacts repositories create "sv"
  --repository-format "docker" \
  --location "australia-southeast1" \
  --description "SV Dockers" \
  --project "<PROJECT_ID>"
```

- build Docker image and upload to specified AR repo under specified tag
  (assuming the Dockerfile is in the current working directory):

```shell
gcloud builds submit \
  --tag australia-southeast1-docker.pkg.dev/PROJECT_ID/sv/IMG:TAG
```
