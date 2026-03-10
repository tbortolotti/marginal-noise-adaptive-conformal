#!/bin/bash

# Slurm parameters
MEMO=20G                             # Memory required (20 GB)
TIME=00-03:00:00                    # Time required (1 h 20 m)
CORE=1                              # Cores required (1)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/exp1000"

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/exp1000"

JOBN="exp1000/prova"
OUT_FILE=$OUT_DIR"/"$JOBN".txt"
COMPLETE=0
if [[ -f $OUT_FILE ]]; then
  COMPLETE=1
fi

if [[ $COMPLETE -eq 0 ]]; then
  # Script to be run
  SCRIPT="exp_cifar_get_features.sh"
  # Define job name
  OUTF=$LOGS"/"$JOBN".out"
  ERRF=$LOGS"/"$JOBN".err"
  # Assemble slurm order for this job
  ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
  # Print order
  echo $ORD
  # Submit order
  $ORD
  # Run command now
  #./$SCRIPT
fi
