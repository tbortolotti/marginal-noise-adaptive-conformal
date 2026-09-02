#!/bin/bash
# Collect the per-repetition experiment outputs produced on a remote cluster into
# the local experiments/results_hpc/ tree, which experiments/make_plots.R reads.
#
# Edit REMOTE below for your setup (host + path to experiments/results/ on the
# machine where the jobs ran).

set -eu

REMOTE="user@cluster.example.edu:/path/to/marginal-noise-adaptive-conformal/experiments/results/"

mkdir -p results_hpc
rsync -auv "$REMOTE"/* results_hpc/
