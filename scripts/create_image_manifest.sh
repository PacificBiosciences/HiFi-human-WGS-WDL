#!/usr/bin/env bash
set -e

# This script generates a manifest file that lists all the docker images used in the WDL files.

find . -type f -name "*.wdl" -exec grep '@sha' -h -r {} \; \
| tr --squeeze-repeats ' ' \
| cut --fields=3 --delimiter=' ' \
| sed 's!~{runtime_attributes.container_registry}!quay.io/pacbio!;s/"//g;' \
| sort --unique \
> ./image_manifest.txt

deepvariant_version=1.10.0
echo "google/deepvariant:${deepvariant_version}" >> ./image_manifest.txt
echo "google/deepvariant:${deepvariant_version}-gpu" >> ./image_manifest.txt

deepsomatic_version=1.10.0
echo "google/deepsomatic:${deepsomatic_version}" >> ./image_manifest.txt
echo "google/deepsomatic:${deepsomatic_version}-gpu" >> ./image_manifest.txt

parabricks=4.7.0-1
echo "nvcr.io/nvidia/clara/clara-parabricks:${parabricks}" >> ./image_manifest.txt

pharmcat_version=2.15.4
echo "pgkb/pharmcat:${pharmcat_version}" >> ./image_manifest.txt
