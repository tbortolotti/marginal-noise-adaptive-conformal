DATA=$1

mkdir -p results_hpc

rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/label_noise_marginal/code_github/experiments/results/* results_hpc/
