#!/bin/bash

#SBATCH --job-name=download_bigearthnet   # Job name
#SBATCH --output=download_bigearthnet.log # Output log file
#SBATCH --error=download_bigearthnet.err  # Error log file
#SBATCH --partition=cpu                   # Use the CPU partition
#SBATCH --cpus-per-task=4                 # Number of CPU cores
#SBATCH --time=05:00:00                   # Max runtime (hh:mm:ss)
#SBATCH --mem=8G                          # Memory allocation

module purge
eval "$(conda shell.bash hook)"
conda activate bigearth

export PYTHONPATH=$(pwd)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Define the Python script to run
cat << EOF > download_dataset.py
import pathlib
from datamodules.bigearthnet_datamodule import download_data

# Define dataset directory and name
dataset_dir = "../datasets"  # Adjust this to your desired path
dataset_name = "bigearthnet-full"

# Download the dataset
download_data(dataset_dir, dataset_name)
EOF

# Run the script
python download_dataset.py
