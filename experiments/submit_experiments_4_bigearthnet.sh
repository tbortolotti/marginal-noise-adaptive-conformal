#!/bin/bash

# Parameters
CONF=201

if [[ $CONF == 201 ]]; then
  # Figure 6
  BATCH_SIZE_LIST=(700 800 900 1000 2000 3000 5000 10000 15000 20000)
  #BATCH_SIZE_LIST=(700)
  ESTIMATE_LIST=("none")
  SEED_LIST=$(seq 1 30)
  #SEED_LIST=$(seq 1 1)

fi

# Slurm parameters
MEMO=64G                             # Memory required (32 GB)
TIME=00-14:00:00                    # Time required (5 h 00 m)
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
  for BATCH_SIZE in "${BATCH_SIZE_LIST[@]}"; do
    for ESTIMATE in "${ESTIMATE_LIST[@]}"; do
			  JOBN="exp"$CONF"/bigearthnet_n"$BATCH_SIZE
			  JOBN=$JOBN"_est"$ESTIMATE_LIST"_"$SEED
			  OUT_FILE=$OUT_DIR"/"$JOBN".txt"
			  COMPLETE=0
			  if [[ -f $OUT_FILE ]]; then
			      COMPLETE=1
			  fi

			  if [[ $COMPLETE -eq 0 ]]; then
			      # Script to be run
			      SCRIPT="exp_bigearthnet.sh $CONF $BATCH_SIZE $ESTIMATE_LIST $SEED"
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
