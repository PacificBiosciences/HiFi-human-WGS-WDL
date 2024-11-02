#!/usr/bin/env bash
set -e

# This script generates a manifest file that lists all the docker images used in the WDL files.

grep '@sha' -h -r workflows/ \
| tr --squeeze-repeats ' ' \
| cut --fields=3 --delimiter=' ' \
| sed 's!~{runtime_attributes.container_registry}!quay.io/pacbio!;s/"//g;' \
| sort --unique \
> ./image_manifest.txt

deepvariant_version=$(grep -m1 'String deepvariant_version' workflows/singleton.wdl | tr -s ' ' | cut -f5 -d' ' | sed 's/"//g')
echo "google/deepvariant:${deepvariant_version}" >> ./image_manifest.txt
echo "google/deepvariant:${deepvariant_version}-gpu" >> ./image_manifest.txt

pharmcat_version=$(grep -m1 'String pharmcat_version' workflows/singleton.wdl | tr -s ' ' | cut -f5 -d' ' | sed 's/"//g')
echo "pgkb/pharmcat:${pharmcat_version}" >> ./image_manifest.txt