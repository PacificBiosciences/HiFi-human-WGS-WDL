#!/usr/bin/env bash
set -e

# create a README.md file from the Home.md file in the wiki directory
# with correct relative links
sed 's!(\./docs/!(./!g;s!\.\./\.\./\.\.!../..!g;' README.md > "$1"/Home.md
cp docs/*.md "$1"/