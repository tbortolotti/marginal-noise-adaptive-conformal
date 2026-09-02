#!/bin/bash
# Slurm submission wrapper for the BigEarthNet download (> 30 GB).
# The simple entry point is  data/download_bigearthnet.sh  from the repo root;
# use this one to run the download as a batch job on a cluster.

#SBATCH --job-name=download_bigearthnet
#SBATCH --output=download_bigearthnet.log
#SBATCH --error=download_bigearthnet.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --time=05:00:00
#SBATCH --mem=8G

set -eu

module purge 2>/dev/null || true
eval "$(conda shell.bash hook)"
conda activate bigearth

export PYTHONPATH="$(pwd)"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# download_dataset.py puts the data in the top-level data/ folder.
python download_dataset.py
