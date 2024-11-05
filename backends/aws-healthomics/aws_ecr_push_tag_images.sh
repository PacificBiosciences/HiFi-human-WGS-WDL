#!/bin/bash

# Author: Heather Ward <heather@dnastack.com>
# License: GPLv2
#
# This script is licensed under the GNU General Public License v2.0.
# See https://www.gnu.org/licenses/old-licenses/gpl-2.0.html for details.

set -euo pipefail

usage() {
cat << EOF

    Usage: $0 -d docker_manifest -r aws_region -a aws_account [OPTIONS]
    Ensure your AWS_PROFILE is set to the correct AWS account (the one you're pushing images to)

    OPTIONS
        -h  Display this message and exit
        -d  Docker manifest; images to pull, tag, and push to the specified ECR, one per line
        -r  AWS region to host images in
        -a  AWS account where images are hosted

EOF
}

log() {
  MESSAGE=$1
  TYPE=${2:-info}

  if [[ "$TYPE" == "err" ]]; then
    >&2 echo "${MESSAGE}"
  elif [[ "$TYPE" == "info" ]]; then
    echo "${MESSAGE}"
  fi
}

push_aws() {
    TARGET_IMAGE=$1
    POLICY_JSON=$2
    AWS_REGION=$3

    repository_name=$(echo "${TARGET_IMAGE}" | cut -d / -f 2- | cut -d ":" -f 1 | cut -d "@" -f 1)
    # Create the repositroy if it does not exist
    if ! aws ecr describe-repositories --repository-names "${repository_name}" --region "${AWS_REGION}" 2> /dev/null; then
        log "Creating repository ${repository_name}"
        aws ecr create-repository \
            --repository-name "${repository_name}" \
            --region "${AWS_REGION}" \
        2> /dev/null
    fi
    log "Setting ECR cross-account policy on ${repository_name}"
    aws ecr set-repository-policy \
        --repository-name "${repository_name}" \
        --policy-text "${POLICY_JSON}" \
        --region "${AWS_REGION}" \
    > /dev/null

    log "Pushing image ${TARGET_IMAGE}"
    docker push "${TARGET_IMAGE}"
}

while getopts "hd:r:a:" OPTION; do
    case "${OPTION}" in
        h) usage; exit;;
        d) DOCKER_MANIFEST="${OPTARG}";;
        r) AWS_REGION="${OPTARG}";;
        a) AWS_ACCOUNT="${OPTARG}";;
        \?) usage; exit;;
    esac
done

DOCKER_MANIFEST=${DOCKER_MANIFEST:-}
AWS_REGION=${AWS_REGION:-}
AWS_ACCOUNT=${AWS_ACCOUNT:-}
POLICY_JSON='{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "OmicsAccessPrincipal",
            "Effect": "Allow",
            "Principal": {
                "Service": "omics.amazonaws.com"
            },
            "Action": [
                "ecr:BatchCheckLayerAvailability",
                "ecr:BatchGetImage",
                "ecr:GetDownloadUrlForLayer"
            ]
        }
    ]
}'

if [[ -z "${DOCKER_MANIFEST}" ]]; then
    usage
    log "Must set -d docker_manifest" err
    exit 1
fi

if [[ -z "${AWS_REGION}" ]]; then
    usage
    log "Must set -r aws_region" err
    exit 1
fi

# Pull, rename, push, and apply policy to all images in the manifest
AWS_CONTAINER_REGISTRY="${AWS_ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com"
aws ecr get-login-password --region "${AWS_REGION}" | docker login --username AWS --password-stdin "${AWS_CONTAINER_REGISTRY}"

while read -r image || [[ -n "${image}" ]]; do
    log "Pulling image ${image}"
    docker pull "${image}"
    # If the image is a @sha256 rather than a tag, tag the image with the sha hash (removing @sha256)
    target_image=${AWS_CONTAINER_REGISTRY}/$(echo "${image}" | awk -F "/" '{print $NF}' | sed "s~@sha256~~")
    log "${target_image}"
    log "Tagging image ${image} as ${target_image}"
    docker tag "${image}" "${target_image}"
    push_aws "${target_image}" "${POLICY_JSON}" "${AWS_REGION}"
done < "${DOCKER_MANIFEST}"
