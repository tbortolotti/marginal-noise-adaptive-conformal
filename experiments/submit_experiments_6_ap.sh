#!/bin/bash

# Parameters
CONF=604

if [[ $CONF == 600 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic1")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  SIGNAL_LIST=(1.0)
  PI_LIST=(1.0)
  C_SCALE_LIST=(1)
  FLIPY_LIST=(0)
  EPSILON_LIST=(0.1)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(10000)
  N_TRAIN2_LIST=(10000)
  N_CAL_LIST=(2000)
  SEED_LIST=(1)

elif [[ $CONF == 601 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic1_easy")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  SIGNAL_LIST=(0.7 1.0 2.0)
  PI_LIST=(1.0)
  C_SCALE_LIST=(1)
  FLIPY_LIST=(0 0.01)
  EPSILON_LIST=(0.1)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(10000)
  N_TRAIN2_LIST=(10000)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 602 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic1")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  SIGNAL_LIST=(0.7 1.0 2.0)
  PI_LIST=(1.0)
  C_SCALE_LIST=(1)
  FLIPY_LIST=(0 0.01)
  EPSILON_LIST=(0.1)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(5000)
  N_TRAIN2_LIST=(5000)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 603 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic1")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  SIGNAL_LIST=(1.0)
  PI_LIST=(1.0)
  C_SCALE_LIST=(1)
  FLIPY_LIST=(0 0.01)
  EPSILON_LIST=(0.1)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(10000)
  N_TRAIN2_LIST=(100 500 1000 5000 10000)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 604 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("syntheticAP")
  NUM_VAR_LIST=(2)
  K_LIST=(4)
  SIGNAL_LIST=(1.0)
  PI_LIST=(1.0)
  C_SCALE_LIST=(1)
  FLIPY_LIST=(0)
  EPSILON_LIST=(0.1)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(500 1000 5000 10000)
  N_TRAIN2_LIST=(500 1000 5000 10000)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 605 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("syntheticAP")
  NUM_VAR_LIST=(2)
  K_LIST=(4)
  SIGNAL_LIST=(1.0)
  PI_LIST=(1.0)
  C_SCALE_LIST=(0.5 0.75 1)
  FLIPY_LIST=(0)
  EPSILON_LIST=(0.1 0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(5000)
  N_TRAIN2_LIST=(5000)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

fi


# Slurm parameters
MEMO=4G                             # Memory required (1 GB)
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
  for MODEL in "${MODEL_LIST[@]}"; do
    for DATA in "${DATA_LIST[@]}"; do
      for NUM_VAR in "${NUM_VAR_LIST[@]}"; do
        for K in "${K_LIST[@]}"; do
          for SIGNAL in "${SIGNAL_LIST[@]}"; do
            for PI in "${PI_LIST[@]}"; do
              for C_SCALE in "${C_SCALE_LIST[@]}"; do
                for FLIPY in "${FLIPY_LIST[@]}"; do
                  for EPSILON in "${EPSILON_LIST[@]}"; do
                    for CONTAMINATION in "${CONTAMINATION_LIST[@]}"; do
                      for N_TRAIN1 in "${N_TRAIN1_LIST[@]}"; do
                        for N_TRAIN2 in "${N_TRAIN2_LIST[@]}"; do
                          for N_CAL in "${N_CAL_LIST[@]}"; do
                            JOBN="exp"$CONF"/"$DATA"_p"$NUM_VAR"_K"$K"_signal"$SIGNAL"_PI"$PI"_CSCALE"$C_SCALE"_"$MODEL"_flipy"$FLIPY"_eps"$EPSILON"_"$CONTAMINATION"_nt1_"$N_TRAIN1"_nt2_"$N_TRAIN2"_nc"$N_CAL"_seed"$SEED
                            OUT_FILE=$OUT_DIR"/"$JOBN".txt"
                            COMPLETE=0
                            #  ls $OUT_FILE
                            if [[ -f $OUT_FILE ]]; then
                                  COMPLETE=1
                            fi

                            if [[ $COMPLETE -eq 0 ]]; then
                              # Script to be run
                              SCRIPT="exp_ap.sh $CONF $MODEL $DATA $NUM_VAR $K $SIGNAL $PI $C_SCALE $FLIPY $EPSILON $CONTAMINATION $N_TRAIN1 $N_TRAIN2 $N_CAL $SEED"
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
            done
          done
        done
      done
    done
  done
done
