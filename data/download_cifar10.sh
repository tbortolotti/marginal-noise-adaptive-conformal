#!/bin/bash
# Download CIFAR-10 into data/cifar10/ (~170 MB) using torchvision.
# Produces  data/cifar10/cifar-10-batches-py/  and  data/cifar10/cifar-10-python.tar.gz
#
# Requires a Python with torch + torchvision on the path (the 'default'
# environment from metadata/dependencies.md).

set -eu
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$HERE/cifar10"
mkdir -p "$DEST"

python - "$DEST" <<'PY'
import sys
from torchvision.datasets import CIFAR10
root = sys.argv[1]
CIFAR10(root=root, train=True, download=True)
CIFAR10(root=root, train=False, download=True)
print("CIFAR-10 ready at", root)
PY
