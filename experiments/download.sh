DATA=$1

mkdir -p results_hpc

#rsync -auv sesia@discovery.usc.edu:/home1/sesia/Workspace/label-noise/code/experiments/results/* results_hpc/

rsync -auv tb_214@discovery.usc.edu:/home1/tb_214/code/marginal-noise-adaptive-conformal/experiments/results/* results_hpc/
