#!/bin/bash

# Parameters
CONF=301

if [[ $CONF == 301 ]]; then
  DATA_LIST=("synthetic4")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  SIGNAL_LIST=(1.0)
  MODEL_LIST=('RFC')
  EPSILON_LIST=(0.1)
  NU_LIST=(0.2)
  CONTAMINATION_LIST=("RRB")
  N_TRAIN_LIST=(10000)
  N_CAL_LIST=(500 1000 2000 5000 10000 20000 50000 100000)
  ESTIMATE_LIST=("none")
  IMB_LIST=(0 0.5 1 2)
  SEED_LIST=$(seq 1 5)
fi


# Slurm parameters
MEMO=5G                             # Memory required (1 GB)
TIME=00-01:00:00                    # Time required (20 m)
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
        for SIGNAL in "${SIGNAL_LIST[@]}"; do
          for MODEL in "${MODEL_LIST[@]}"; do
            for EPSILON in "${EPSILON_LIST[@]}"; do
	            for NU in "${NU_LIST[@]}"; do
       		      for CONTAMINATION in "${CONTAMINATION_LIST[@]}"; do
                  for N_TRAIN in "${N_TRAIN_LIST[@]}"; do
                    for N_CAL in "${N_CAL_LIST[@]}"; do
                      for ESTIMATE in "${ESTIMATE_LIST[@]}"; do
                        for IMB in "${IMB_LIST[@]}"; do
                          JOBN="exp"$CONF"/"$DATA"_p"$NUM_VAR"_K"$K"_signal"$SIGNAL"_"$MODEL"_eps"$EPSILON"_nu"$NU"_"$CONTAMINATION"_nt"$N_TRAIN"_nc"$N_CAL"_est"$ESTIMATE"_imb"$IMB"_seed"$SEED
                          OUT_FILE=$OUT_DIR"/"$JOBN".txt"
                          COMPLETE=0
                          #                    ls $OUT_FILE
                          if [[ -f $OUT_FILE ]]; then
                              COMPLETE=1
			                    fi

			                    if [[ $COMPLETE -eq 0 ]]; then
                              # Script to be run
                              SCRIPT="exp_classimbalance.sh $CONF $DATA $NUM_VAR $K $SIGNAL $MODEL $EPSILON $NU $CONTAMINATION $N_TRAIN $N_CAL $ESTIMATE $IMB $SEED"
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
