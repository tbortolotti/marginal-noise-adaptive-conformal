#!/bin/bash

module purge
eval "$(conda shell.bash hook)"
conda activate default

export OPENBLAS_NUM_THREADS=1

python3 exp_cifar_ap.py $1 $2 $3 $4 $5 $6 $7 $8
