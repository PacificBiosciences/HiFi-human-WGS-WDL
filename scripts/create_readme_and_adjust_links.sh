#!/usr/bin/env bash
set -e

# create a README.md file from the Home.md file in the wiki directory
# with correct relative links
sed 's!(\./!(./wiki/!g;s!\.\./\.\.!../../..!g;' ./wiki/Home.md > README.md
