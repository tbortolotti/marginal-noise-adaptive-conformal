#!/bin/bash

# Parameters
CONF=1103

if [[ $CONF == 1100 ]]; then
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("real")
  N_TRAIN_LIST=(5000)
  N_CLEAN_LIST=(500)
  N_CAL_LIST=(1000)
  CONTAMINATION_EXP_FLAG="true"
  SEED_LIST=(1)

elif [[ $CONF == 1101 ]]; then
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN_LIST=(5000)
  N_CLEAN_LIST=(500)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000)
  CONTAMINATION_EXP_FLAG="false"
  SEED_LIST=$(seq 1 20)

elif [[ $CONF == 1102 ]]; then
  EPSILON_LIST=(0.1)
  NU_LIST=(0.2)
  CONTAMINATION_LIST=("uniform" "block" "RRB")
  N_TRAIN_LIST=(5000)
  N_CLEAN_LIST=(500)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000)
  CONTAMINATION_EXP_FLAG="true"
  SEED_LIST=$(seq 1 20)

elif [[ $CONF == 1103 ]]; then
  EPSILON_LIST=(0.016)
  NU_LIST=(0)
  CONTAMINATION_LIST=("real")
  N_TRAIN_LIST=(5000)
  N_CLEAN_LIST=(500)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000)
  CONTAMINATION_EXP_FLAG="true"
  SEED_LIST=$(seq 1 20)

fi


# Slurm parameters
#MEMO=64G
#TIME=00-04:00:00
#CORE=1

# Assemble order prefix
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=main"

MEMO=64G 
TIME=00-06:00:00
CORE=1
#ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME" --account=sesia_1124 --partition=gpu --gres=gpu:p100:1"
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

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
        for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
          for N_CLEAN in "${N_CLEAN_LIST[@]}"; do
            for N_CAL in "${N_CAL_LIST[@]}"; do
              JOBN="exp"$CONF"/bigearthnet_eps"$EPSILON
              JOBN=$JOBN"_nu"$NU"_"$CONTAMINATION"_nt"$N_TRAIN"_ncl"$N_CLEAN"_nc"$N_CAL"_seed"$SEED
              OUT_FILE=$OUT_DIR"/"$JOBN".txt"
              COMPLETE=0
              #  ls $OUT_FILE
              if [[ -f $OUT_FILE ]]; then
                    COMPLETE=1
              fi

              if [[ $COMPLETE -eq 0 ]]; then
                # Script to be run
                SCRIPT="exp_bigearthnet_EM_NN.sh $CONF $EPSILON $NU $CONTAMINATION $N_TRAIN $N_CLEAN $N_CAL $CONTAMINATION_EXP_FLAG $SEED"
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
