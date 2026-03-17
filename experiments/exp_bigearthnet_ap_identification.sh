#!/bin/bash

module purge
eval "$(conda shell.bash hook)"
conda activate default

export OPENBLAS_NUM_THREADS=1

python3 exp_bigearthnet_ap_identification.py $1 $2 $3 $4 $5 $6 $7
