#!/bin/bash
# Download the pre-processed BigEarthNet ("hub" format) archive into
# data/bigearthnet-full/  -- this is the data used by the BigEarthNet experiment
# (experiments/exp_bigearthnet.py).
#
#   ############################################################
#   #  THIS DOWNLOAD IS LARGER THAN 30 GB (compressed archive  #
#   #  ~30 GB, more once extracted).  Run it on a cluster or a #
#   #  machine with plenty of free disk and a stable network   #
#   #  connection -- NOT a laptop.                             #
#   ############################################################
#
# It is NOT needed to rebuild the paper figures: those use the pre-computed
# results shipped in experiments/results/. It is only needed to re-run the
# BigEarthNet experiment from scratch.
#
# Requires `gdown` (pip install gdown; already in the 'bigearth' environment).
# The Google Drive id below is the "bigearthnet-full" entry of GDRIVE_URLS in
# third_party/bigearthnet/datamodules/bigearthnet_datamodule.py.

set -eu
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" != "--yes" ]; then
  echo "BigEarthNet is > 30 GB. This will download into $HERE/bigearthnet-full/."
  echo "Re-run to confirm:   ./download_bigearthnet.sh --yes"
  exit 1
fi

GDRIVE_URL="https://drive.google.com/file/d/1isUcPQvCn1xc5GWEmPDOqj_l24QH2ZtA/view?usp=sharing"
TAR="$HERE/bigearthnet-full.tar"
DEST="$HERE/bigearthnet-full"

if [ -d "$DEST" ]; then
  echo "Already present at $DEST -- nothing to do."
  exit 0
fi

if [ ! -f "$TAR" ]; then
  python - "$GDRIVE_URL" "$TAR" <<'PY'
import sys, gdown
gdown.download(sys.argv[1], sys.argv[2], fuzzy=True)
PY
fi

echo "Extracting $TAR ..."
tar -xf "$TAR" -C "$HERE"
echo "BigEarthNet ready at $DEST/"
echo "(you may now delete $TAR to reclaim space)"
