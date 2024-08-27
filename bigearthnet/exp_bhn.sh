#!/bin/bash
#SBATCH --job-name=resnet50_training
#SBATCH --output=resnet50_training_%j.out  # Output file (%j will be replaced with the job ID)
#SBATCH --error=resnet50_training_%j.err   # Error file (%j will be replaced with the job ID)
#SBATCH --ntasks=1                        # Number of tasks (usually 1 for a single script)
#SBATCH --cpus-per-task=4                 # Number of CPU cores per task
#SBATCH --mem=32G                         # Increase memory allocation to 32 GB
#SBATCH --time=24:00:00                   # Increase time limit to 24 hours

module purge
eval "$(conda shell.bash hook)"
conda activate default

export OPENBLAS_NUM_THREADS=1

python3 BHN_ResNet_train.py
