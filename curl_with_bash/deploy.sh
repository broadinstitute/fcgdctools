#!/bin/bash

# method name
METHOD_NAME=curl
VERSION=7.75.0

# build method docker
docker build -t gcr.io/broad-getzlab-workflows/$METHOD_NAME:${VERSION} .
docker push gcr.io/broad-getzlab-workflows/$METHOD_NAME:${VERSION}
