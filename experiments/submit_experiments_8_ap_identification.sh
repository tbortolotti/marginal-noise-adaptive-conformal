#!/bin/bash

# Parameters
CONF=801

if [[ $CONF == 800 ]]; then
  DATA_LIST=("syntheticAP")
  NUM_VAR_LIST=(2)
  K_LIST=(4)
  PI_LIST=(1)
  DSHIFT_LIST=(0)
  C_SCALE_LIST=(1)
  SEASY_LIST=(1)
  SHARD_LIST=(1)
  REASY_LIST=(1.5)
  RHARD_LIST=(4.0)
  FLIPY_LIST=(0)
  MODEL_LIST=('RFC')
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(10000)
  N_TRAIN2_LIST=(1000)
  SEED_LIST=(1)

elif [[ $CONF == 801 ]]; then
  DATA_LIST=("syntheticAP")
  NUM_VAR_LIST=(2)
  K_LIST=(4)
  PI_LIST=(0 0.25 0.5 0.75 1)
  DSHIFT_LIST=(0)
  C_SCALE_LIST=(1)
  SEASY_LIST=(1)
  SHARD_LIST=(1)
  REASY_LIST=(1.5)
  RHARD_LIST=(4.0)
  FLIPY_LIST=(0)
  MODEL_LIST=('RFC')
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(1000 5000 10000 20000 50000 100000)
  N_TRAIN2_LIST=(10000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 802 ]]; then
  DATA_LIST=("syntheticAP")
  NUM_VAR_LIST=(2)
  K_LIST=(4)
  PI_LIST=(0 0.25 0.5 0.75 1)
  DSHIFT_LIST=(0)
  C_SCALE_LIST=(1)
  SEASY_LIST=(1)
  SHARD_LIST=(1)
  REASY_LIST=(1.5)
  RHARD_LIST=(4.0)
  FLIPY_LIST=(0)
  MODEL_LIST=('RFC')
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN1_LIST=(10000)
  N_TRAIN2_LIST=(500 1000 2000 5000 10000 20000 50000 100000)
  SEED_LIST=$(seq 1 5)
fi


# Slurm parameters
MEMO=5G                             # Memory required (1 GB)
TIME=00-08:00:00                    # Time required (4 h)
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
    for NUM_VAR in "${NUM_VAR_LIST[@]}"; do
      for K in "${K_LIST[@]}"; do
        for PI in "${PI_LIST[@]}"; do
          for DSHIFT in "${DSHIFT_LIST[@]}"; do
            for C_SCALE in "${C_SCALE_LIST[@]}"; do
              for SEASY in "${SEASY_LIST[@]}"; do
                for SHARD in "${SHARD_LIST[@]}"; do
                  for REASY in "${REASY_LIST[@]}"; do
                    for RHARD in "${RHARD_LIST[@]}"; do
                      for FLIPY in "${FLIPY_LIST[@]}"; do
                        for MODEL in "${MODEL_LIST[@]}"; do
                          for EPSILON in "${EPSILON_LIST[@]}"; do
                            for CONTAMINATION in "${CONTAMINATION_LIST[@]}"; do
                              for N_TRAIN1 in "${N_TRAIN1_LIST[@]}"; do
                                for N_TRAIN2 in "${N_TRAIN2_LIST[@]}"; do
                                  JOBN="exp"$CONF"/"$DATA"_p"$NUM_VAR"_K"$K"_PI"$PI"_DSHIFT"$DSHIFT"_CSCALE"$C_SCALE"_SEASY"$SEASY"_SHARD"$SHARD"_REASY"$REASY"_RHARD"$RHARD"_flipy"$FLIPY"_"$MODEL"_eps"$EPSILON"_"$CONTAMINATION"_nt1_"$N_TRAIN1"_nt2_"$N_TRAIN2"_seed"$SEED
                                  OUT_FILE=$OUT_DIR"/"$JOBN".txt"
                                  COMPLETE=0
                                  #  ls $OUT_FILE
                                  if [[ -f $OUT_FILE ]]; then
                                        COMPLETE=1
                                  fi

                                  if [[ $COMPLETE -eq 0 ]]; then
                                    # Script to be run
                                    SCRIPT="exp_ap_identification.sh $CONF $DATA $NUM_VAR $K $PI $DSHIFT $C_SCALE $SEASY $SHARD $REASY $RHARD $FLIPY $MODEL $EPSILON $CONTAMINATION $N_TRAIN1 $N_TRAIN2 $SEED"
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
    done
  done
done
