# Dockerfiles

- [Dockerfiles](#dockerfiles)
  - [Create AR repo](#create-ar-repo)
  - [Build Docker image and upload to AR/GCR repo](#build-docker-image-and-upload-to-argcr-repo)
  - [List AR Docker images](#list-ar-docker-images)
  - [Copy Docker images across GCP](#copy-docker-images-across-gcp)

We use [Cloud Build](https://cloud.google.com/build/docs/quickstart-build) to
build all Docker images for the main structural variant tools; these then get
automatically uploaded to
[Artifact Registry](https://cloud.google.com/artifact-registry/docs/docker/quickstart)
(AR) (or
[Google Container Registry](https://cloud.google.com/container-registry) (GCR)
since Cromwell does not support AR as of v56).

## Create AR repo

- Create AR repo named `sv` under a given Google Cloud project (done once):

```shell
gcloud artifacts repositories create "sv"
  --repository-format "docker" \
  --location "australia-southeast1" \
  --description "SV Dockers" \
  --project "<PROJECT_ID>"
```

## Build Docker image and upload to AR/GCR repo

- Build Docker image and upload to specified AR/GCR repo under specified tag
  (assuming the Dockerfile is in the current working directory):

```shell
gcloud builds submit \
  --tag australia-southeast1-docker.pkg.dev/PROJECT_ID/sv/IMG:TAG
  #--tag us.gcr.io/PROJECT_ID/IMG:TAG
```

## List AR Docker images

- List AR Docker images:

```shell
gcloud artifacts docker images list \
  australia-southeast1-docker.pkg.dev/peter-dev-302805/sv/ \
  --include-tags
```

## Copy Docker images across GCP

- For other Docker images where Dockerfiles are out of sync with tagged
  versions, we simply pull the tagged version and push to our own AR/GCR repo:

```shell
docker pull $TAG
docker tag $TAG $OUR_TAG
docker push $OUR_TAG
```
