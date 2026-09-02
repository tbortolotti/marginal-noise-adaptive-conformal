#!/bin/bash

# On a Slurm cluster this selects the right environment; a no-op elsewhere.
module purge 2>/dev/null || true
if command -v conda >/dev/null 2>&1; then
  eval "$(conda shell.bash hook)"
  conda activate default 2>/dev/null || true   # edit if your env name differs
fi

export OPENBLAS_NUM_THREADS=1

python3 exp_T_estimation.py $1 $2 $3 $4 $5 $6 $7 $8 $9 "${10}" "${11}" "${12}" "${13}"
