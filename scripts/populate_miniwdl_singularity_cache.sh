#!/usr/bin/env bash
set -eo pipefail

USAGE="Given a manifest file with docker images, this script populates the Singularity cache with those images.
Usage: $0 <image_manifest_file> <miniwdl_singularity_cache_dir>"


# Check if the first argument is -h or --help
if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  echo -e "${USAGE}"
  exit 0
fi
# Check if at least two arguments are provided
if [ $# -lt 2 ]; then
  echo -e "${USAGE}"
  exit 1
fi

image_manifest_file=$1
miniwdl_singularity_cache_dir=$2
# Check if the image manifest file exists and is readable, and if the cache directory exists and is writable
[ -r "${image_manifest_file}" ] || (echo "${image_manifest_file} is not readable." >&2 && exit 1)
if [ ! -d "${miniwdl_singularity_cache_dir}" ]; then
  echo "${miniwdl_singularity_cache_dir} does not exist.  Creating it now..."
  mkdir -p "${miniwdl_singularity_cache_dir}" || (echo "Could not create ${miniwdl_singularity_cache_dir}/." >&2 && exit 1)
fi
[ -w "${miniwdl_singularity_cache_dir}" ] || (echo "${miniwdl_singularity_cache_dir}/ is not writable" >&2 && exit 1)
singularity --version || (echo "singularity is not in path. Please install Singularity to use this script." >&2 && exit 1)

# manifest file should contain one image per line, with no empty lines
# image lines should be in the format: <image_name>:<tag> or <image_name>@sha256:<digest>
# e.g., google/deepvariant:1.9.0 or quay.io/pacbio/some_image@sha256:abc123...
while read -r image; do
  if [[ -n "${image}" ]]; then
    image_url="docker://${image}"
    # miniwdl singularity backend replaces ':' and '/' with '_' in the SIF file name
    sif_path="${image_url//:/_}"
    sif_path="${sif_path//\//_}"
    sif_path="${miniwdl_singularity_cache_dir}/${sif_path}.sif"
    if [ -f "${sif_path}" ]; then
      echo "Singularity image already exists: ${sif_path}" >&2
    else
      echo "Pulling Singularity image: ${image_url}"
      singularity pull "${sif_path}" "${image_url}"
    fi
  fi
done < "${image_manifest_file}"
