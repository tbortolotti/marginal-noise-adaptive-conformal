#!/bin/bash

DATA=$1

# Get the local system's username
local_username=$(whoami)

mkdir -p results_hpc

# Check the local username and set the remote username accordingly
if [ "$local_username" = "msesia" ]; then
    rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/label_noise_marginal/code_github/experiments/results/* results_hpc/
elif [ "$local_username" = "tb_214" ]; then
    rsync -auv tb_214@discovery.usc.edu:/home1/tb_214/code/marginal-noise-adaptive-conformal/experiments/results/* results_hpc/
elif [ "$local_username" = "teresabortolotti" ]; then
    rsync -auv tb_214@discovery.usc.edu:/home1/tb_214/code/marginal-noise-adaptive-conformal/experiments/results/* results_hpc/
else
    # Print an error message if the user is not recognized
    echo "Error: Unknown user '$local_username'. We don't know what to do for this user."
    exit 1
fi
