#!/bin/bash
# Download the CIFAR-10H human-annotation counts into data/cifar10h/ (~800 KB).
# Produces  data/cifar10h/cifar10h-counts.npy
#
# Source: https://github.com/jcpeterson/cifar-10h  (Peterson et al., 2019)

set -eu
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$HERE/cifar10h"
mkdir -p "$DEST"

URL="https://raw.githubusercontent.com/jcpeterson/cifar-10h/master/data/cifar10h-counts.npy"
curl -fL "$URL" -o "$DEST/cifar10h-counts.npy"
echo "CIFAR-10H ready at $DEST/cifar10h-counts.npy"
