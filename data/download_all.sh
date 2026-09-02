#!/bin/bash
# Download the real datasets used by the paper into this data/ folder.
#
#   CIFAR-10    -> data/cifar10/            (~170 MB)   torchvision
#   CIFAR-10H   -> data/cifar10h/           (~1 MB)     github.com/jcpeterson/cifar-10h
#   BigEarthNet -> data/bigearthnet-full/   (> 30 GB)   run with --with-bigearthnet
#
# The synthetic-data experiments need none of this. The figure-rebuilding step
# (Step 2 of the workflow) needs none of this either -- it uses the shipped
# results. Only re-running the CIFAR-10 / BigEarthNet experiments needs the data.

set -eu
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$HERE/download_cifar10.sh"
"$HERE/download_cifar10h.sh"

if [ "${1:-}" = "--with-bigearthnet" ]; then
  "$HERE/download_bigearthnet.sh" --yes
else
  echo
  echo "Skipped BigEarthNet (> 30 GB). To fetch it:"
  echo "    ./download_bigearthnet.sh --yes"
  echo "  or re-run this script as:"
  echo "    ./download_all.sh --with-bigearthnet"
fi
