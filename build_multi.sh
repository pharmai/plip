#!/bin/bash

version=$(grep version plip/basic/config.py | awk '{print $3}' | sed s/\'//g)
echo "Building multi-architecture version $version..."
export DOCKER_CLI_EXPERIMENTAL=enabled
docker run --rm --privileged docker/binfmt:820fdd95a9972a5308930a2bdfb8573dd4447ad3
docker buildx create --use --name mybuilder
docker buildx build -t pharmai/plip:"$version"-multi --platform=linux/arm/v7,linux/arm64/v8,linux/amd64 --push .