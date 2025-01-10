#!/bin/bash

module purge
eval "$(conda shell.bash hook)"
conda activate bigearth

export PYTHONPATH=$(pwd)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

python exp_bigearthnet_lc.py $1 $2 $3 $4 $5
