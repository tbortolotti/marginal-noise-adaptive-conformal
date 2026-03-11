#!/bin/bash

# Parameters
CONF=901

if [[ $CONF == 900 ]]; then
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(1000)
  N_TRAIN2_LIST=(1000)
  SEED_LIST=(1)

elif [[ $CONF == 901 ]]; then
  EPSILON_LIST=(0 0.05 0.1 0.2)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(500 1000 2000 5000)
  N_TRAIN2_LIST=(1000)
  SEED_LIST=$(seq 1 10)

elif [[ $CONF == 902 ]]; then
  EPSILON_LIST=(0 0.05 0.1 0.2)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(1000)
  N_TRAIN2_LIST=(500 1000 2000 5000)
  SEED_LIST=$(seq 1 10)

fi


# Slurm parameters
MEMO=64G
TIME=00-03:00:00
CORE=1

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS
mkdir -p $LOGS"/exp"$CONF

OUT_DIR="results"
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR"/exp"$CONF


# Loop over configurations
for SEED in $SEED_LIST; do
  for EPSILON in "${EPSILON_LIST[@]}"; do
    for NU in "${NU_LIST[@]}"; do
      for CONTAMINATION in "${CONTAMINATION_LIST[@]}"; do
        for N_TRAIN1 in "${N_TRAIN1_LIST[@]}"; do
          for N_TRAIN2 in "${N_TRAIN2_LIST[@]}"; do
            JOBN="exp"$CONF"/cifar10_eps"$EPSILON
            JOBN=$JOBN"_nu"$NU"_"$CONTAMINATION
            JOBN=$JOBN"_nt1_"$N_TRAIN1"_nt2_"$N_TRAIN2"_"$SEED
            OUT_FILE=$OUT_DIR"/"$JOBN".txt"
            COMPLETE=0
            if [[ -f $OUT_FILE ]]; then
              COMPLETE=1
            fi

            if [[ $COMPLETE -eq 0 ]]; then
              # Script to be run
              SCRIPT="exp_cifar_ap_identification.sh $CONF $EPSILON $NU $CONTAMINATION $N_TRAIN1 $N_TRAIN2 $SEED"
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
          done
        done
      done
    done
  done
done
