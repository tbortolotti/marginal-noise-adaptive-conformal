#!/bin/bash

# Parameters
CONF=1000

if [[ $CONF == 1000 ]]; then
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(2000)
  N_TRAIN2_LIST=(5000)
  N_CAL_LIST=(10000)
  SEED_LIST=(1)

elif [[ $CONF == 1001 ]]; then
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(2000)
  N_TRAIN2_LIST=(1000 2000 5000)
  N_CAL_LIST=(500 1000 2000 5000 10000)
  SEED_LIST=$(seq 1 7)

fi


# Slurm parameters
#MEMO=64G
#TIME=00-04:00:00
#CORE=1

# Assemble order prefix
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

MEMO=32G 
TIME=00-04:00:00
CORE=1
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=gpu --gres=gpu:p100:1"

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
            for N_CAL in "${N_CAL_LIST[@]}"; do
              JOBN="exp"$CONF"/cifar10_eps"$EPSILON
              JOBN=$JOBN"_nu"$NU"_"$CONTAMINATION"_nt1_"$N_TRAIN1"_nt2_"$N_TRAIN2"_nc"$N_CAL"_seed"$SEED
              OUT_FILE=$OUT_DIR"/"$JOBN".txt"
              COMPLETE=0
              #  ls $OUT_FILE
              if [[ -f $OUT_FILE ]]; then
                    COMPLETE=1
              fi

              if [[ $COMPLETE -eq 0 ]]; then
                # Script to be run
                SCRIPT="exp_cifar_ap.sh $CONF $EPSILON $NU $CONTAMINATION $N_TRAIN1 $N_TRAIN2 $N_CAL $SEED"
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
done
