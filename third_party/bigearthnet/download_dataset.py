"""Download the pre-processed BigEarthNet archive.

The supported entry point is  data/download_bigearthnet.sh  (from the repo root).
This script is the underlying library call; run it from this directory with the
`bigearth` environment active. Both put the data in the top-level data/ folder.
"""
import os
from datamodules.bigearthnet_datamodule import download_data

# Top-level data/ folder, relative to third_party/bigearthnet/
dataset_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")
dataset_name = "bigearthnet-full"

# Download the dataset (> 30 GB)
download_data(dataset_dir, dataset_name)
