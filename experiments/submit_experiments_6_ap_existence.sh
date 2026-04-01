#!/bin/bash

# Parameters
CONF=610

if [[ $CONF == 610 ]]; then
  DATA_LIST=("syntheticAP")
  SCENARIO_LIST=("scenario1")
  EPSILON_LIST=(0.1)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(10000)
  N_TRAIN2_LIST=(1000)
  SEED_LIST=(1)

elif [[ $CONF == 611 ]]; then
  DATA_LIST=("syntheticAP")
  SCENARIO_LIST=("scenario1" "scenario2" "scenario3")
  EPSILON_LIST=(0.1)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(10000)
  N_TRAIN2_LIST=(500 1000 5000 10000)
  SEED_LIST=$(seq 1 5)
fi


# Slurm parameters
MEMO=5G                             # Memory required (1 GB)
TIME=00-04:00:00                    # Time required (4 h)
CORE=1                              # Cores required (1)

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
  for DATA in "${DATA_LIST[@]}"; do
    for SCENARIO in "${SCENARIO_LIST[@]}"; do
      for EPSILON in "${EPSILON_LIST[@]}"; do
        for CONTAMINATION in "${CONTAMINATION_LIST[@]}"; do
          for N_TRAIN1 in "${N_TRAIN1_LIST[@]}"; do
            for N_TRAIN2 in "${N_TRAIN2_LIST[@]}"; do
              JOBN="exp"$CONF"/"$DATA"_"$SCENARIO"_eps"$EPSILON"_"$CONTAMINATION"_nt1_"$N_TRAIN1"_nt2_"$N_TRAIN2"_seed"$SEED
              OUT_FILE=$OUT_DIR"/"$JOBN".txt"
              COMPLETE=0
              #  ls $OUT_FILE
              if [[ -f $OUT_FILE ]]; then
                    COMPLETE=1
              fi

              if [[ $COMPLETE -eq 0 ]]; then
                # Script to be run
                SCRIPT="exp_ap_existence.sh $CONF $DATA $SCENARIO $EPSILON $CONTAMINATION $N_TRAIN1 $N_TRAIN2 $SEED"
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
