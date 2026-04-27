#!/bin/bash

# Parameters
CONF=621

if [[ $CONF == 620 ]]; then
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  N_LIST=(1000)
  N_CLEAN_LIST=(100)
  PI_CLEAN_LIST=(0)
  RANDOM_FLAG_LIST=("true")
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  SEED_LIST=(1)

elif [[ $CONF == 621 ]]; then
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  N_LIST=(500 1000 2000 5000 10000 20000)
  N_CLEAN_LIST=(100 500 1000 5000 10000 20000)
  PI_CLEAN_LIST=(0)
  RANDOM_FLAG_LIST=("false")
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 622 ]]; then
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  N_LIST=(500 1000 2000 5000 10000 20000)
  N_CLEAN_LIST=(0)
  PI_CLEAN_LIST=(0.1 0.2 0.5)
  RANDOM_FLAG_LIST=("false")
  EPSILON_LIST=(0.1)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 623 ]]; then
  DATA_LIST=("synthetic6")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  N_LIST=(500 1000 2000 5000 10000 20000)
  N_CLEAN_LIST=(100)
  PI_CLEAN_LIST=(0)
  RANDOM_FLAG_LIST=("false")
  EPSILON_LIST=(0 0.05 0.1 0.2)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
  SEED_LIST=$(seq 1 5)

elif [[ $CONF == 624 ]]; then
  DATA_LIST=("synthetic1" "synthetic2" "synthetic3")
  NUM_VAR_LIST=(20)
  K_LIST=(4)
  N_LIST=(500 1000 2000 5000 10000 20000 50000)
  N_CLEAN_LIST=(100)
  PI_CLEAN_LIST=(0)
  RANDOM_FLAG_LIST=("false")
  EPSILON_LIST=(0.2)
  NU_LIST=(0)
  CONTAMINATION_LIST=("uniform")
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
    for NUM_VAR in "${NUM_VAR_LIST[@]}"; do
      for K in "${K_LIST[@]}"; do
        for N in "${N_LIST[@]}"; do
          for N_CLEAN in "${N_CLEAN_LIST[@]}"; do
            for PI_CLEAN in "${PI_CLEAN_LIST[@]}"; do
              for RANDOM_FLAG in "${RANDOM_FLAG_LIST[@]}"; do
                for EPSILON in "${EPSILON_LIST[@]}"; do
                  for NU in "${NU_LIST[@]}"; do
                    for CONTAMINATION in "${CONTAMINATION_LIST[@]}"; do
                      JOBN="exp"$CONF"/"$DATA"_p"$NUM_VAR"_K"$K"_n"$N"_ncl"$N_CLEAN"_pic"$PI_CLEAN"_eps"$EPSILON"_nu"$NU"_"$CONTAMINATION"_seed"$SEED
                      OUT_FILE=$OUT_DIR"/"$JOBN".txt"
                      COMPLETE=0
                      #  ls $OUT_FILE
                      if [[ -f $OUT_FILE ]]; then
                            COMPLETE=1
                      fi

                      if [[ $COMPLETE -eq 0 ]]; then
                        # Script to be run
                        SCRIPT="exp_T_estimation_EM_NN.sh $CONF $DATA $NUM_VAR $K $N $N_CLEAN $PI_CLEAN $RANDOM_FLAG $EPSILON $NU $CONTAMINATION $SEED"
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
