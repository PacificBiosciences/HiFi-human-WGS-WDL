#!/usr/bin/env bash
set -e

# create a README.md file from the Home.md file in the wiki directory
# with correct relative links
sed 's!(\./!(./wiki/!g;s!\.\./\.\.!../../..!g;s!\.\./images!./images!g' ./wiki/Home.md > README.md
