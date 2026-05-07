#!/bin/bash

# Parameters
CONF=710

if [[ $CONF == 710 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN_LIST=(1000)
  N_CLEAN_LIST=(100)
  PI_CLEAN_LIST=(0)
  N_CAL_LIST=(2000)
  SEED_LIST=(1)

elif [[ $CONF == 711 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN_LIST=(10000)
  N_CLEAN_LIST=(100 500 1000)
  PI_CLEAN_LIST=(0)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 712 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN_LIST=(10000)
  N_CLEAN_LIST=(0)
  PI_CLEAN_LIST=(0.05 0.1 0.2)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 713 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN_LIST=(1000 5000 10000)
  N_CLEAN_LIST=(100 500)
  PI_CLEAN_LIST=(0)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 714 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN_LIST=(1000 5000 10000)
  N_CLEAN_LIST=(0)
  PI_CLEAN_LIST=(0.1)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000)
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 715 ]]; then
  MODEL_LIST=('RFC')
  DATA_LIST=("synthetic1" "synthetic2" "synthetic3")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  EPSILON_LIST=(0.2)
  CONTAMINATION_LIST=("uniform")
  N_TRAIN_LIST=(10000)
  N_CLEAN_LIST=(500)
  PI_CLEAN_LIST=(0)
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
          for EPSILON in "${EPSILON_LIST[@]}"; do
            for CONTAMINATION in "${CONTAMINATION_LIST[@]}"; do
              for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
                for N_CLEAN in "${N_CLEAN_LIST[@]}"; do
                  for PI_CLEAN in "${PI_CLEAN_LIST[@]}"; do
                    for N_CAL in "${N_CAL_LIST[@]}"; do
                      JOBN="exp"$CONF"/"$DATA"_p"$NUM_VAR"_K"$K"_"$MODEL"_eps"$EPSILON"_"$CONTAMINATION"_nt_"$N_TRAIN"_ncl_"$N_CLEAN"_picl_"$PI_CLEAN"_nc"$N_CAL"_seed"$SEED
                      OUT_FILE=$OUT_DIR"/"$JOBN".txt"
                      COMPLETE=0
                      #  ls $OUT_FILE
                      if [[ -f $OUT_FILE ]]; then
                            COMPLETE=1
                      fi

                      if [[ $COMPLETE -eq 0 ]]; then
                        # Script to be run
                        SCRIPT="exp_EM_NN.sh $CONF $MODEL $DATA $NUM_VAR $K $EPSILON $CONTAMINATION $N_TRAIN $N_CLEAN $PI_CLEAN $N_CAL $SEED"
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
