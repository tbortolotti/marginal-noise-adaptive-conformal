#!/bin/bash

module purge
eval "$(conda shell.bash hook)"
conda activate bigearth

export OPENBLAS_NUM_THREADS=1

python3 exp_bigearthnet_EM_NN.py $1 $2 $3 $4 $5 $6 $7 $8 $9
