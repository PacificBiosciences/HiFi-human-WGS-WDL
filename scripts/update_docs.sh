#!/usr/bin/env bash
set -e

# create a README.md file from the Home.md file in the wiki directory
# with correct relative links
sed 's!(\./!(./docs/!g;s!\.\./\.\.!../../..!g;' "$1"/Home.md > README.md
cp "$1"/*.md docs/
rm docs/Home.md docs/_Sidebar.md