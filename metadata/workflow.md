<!-- Reproducibility workflow. Referenced by README.md. -->

Reproducing a figure involves two stages: **(1)** running the experiments, whose per-
repetition outputs are collected under `experiments/results_hpc/`, and **(2)** building the
figures from those results with `experiments/make_plots.R`.

The experiment outputs (`experiments/results/`, `experiments/results_hpc/`) are **not**
tracked in this repository, because of their size. Stage 1 regenerates them (Slurm cluster).
The pre-computed results also accompany the paper's supplementary material; with them in
place, Stage 2 alone takes a few minutes on a laptop.

## Step 0 -- Install dependencies

See `metadata/dependencies.md`. In short:

```bash
# Python: synthetic-data + CIFAR-10 experiments, the cln package, the notebooks
conda create -n default python=3.11 && conda activate default
pip install -r requirements.txt

# Python: BigEarthNet experiment only (older stack)
conda create -n bigearth python=3.8 && conda activate bigearth
pip install -r requirements-bigearthnet.txt

# R: figures
Rscript -e 'install.packages(c("tidyverse","latex2exp","RColorBrewer","cowplot","patchwork","ggforce"))'
```

The `cln` package is imported by path (`sys.path.append("..")` in the experiment scripts);
no `pip install` of `cln` itself is required. Run all experiment scripts from inside
`experiments/`.

## Step 0b -- Download the real data (only for the CIFAR-10 / BigEarthNet experiments)

Not needed for the synthetic experiments, nor for Stage 2 if you already have the results.
From the repo root:

```bash
cd data
./download_all.sh                     # CIFAR-10 (~170 MB) + CIFAR-10H (~1 MB)
./download_all.sh --with-bigearthnet  # also BigEarthNet
```

**BigEarthNet is larger than 30 GB** -- run `./download_bigearthnet.sh --yes` on a cluster
or a machine with ample free disk, not a laptop. Details and per-dataset scripts:
`data/README.md`. The experiment scripts read from `data/` by default.

## Step 1 -- Run the experiments (Slurm cluster)

**To submit every experiment that backs a manuscript figure in one go:**

```bash
cd experiments
./submit_paper_experiments.sh            # submit all figure experiments
./submit_paper_experiments.sh --dry-run  # print what would run
./submit_paper_experiments.sh --list     # figure -> (launcher, CONF) map
```

This excludes the two Jupyter-notebook studies (appendix figures A1-A6 and A13-A14 -- run
those from `exp_simplified_methods.ipynb` and `exp_richardson.ipynb`). It submits `sbatch`
jobs, so run it from a Slurm head node; completed jobs are skipped, so it is safe to re-run.

**To run one family manually:** each experiment family has a launcher in `experiments/`:

| Launcher | Experiments |
|---|---|
| `submit_experiments_1.sh` | synthetic data: contamination strength / model / number of classes / marginal-vs-conditional |
| `submit_experiments_2_comparisons.sh` | comparison with label-conditional and Clarkson et al. (2024) |
| `submit_experiments_3_T_estimation.sh` | estimation of the noise-transition matrix T |
| `submit_experiments_4.sh` | adaptive method using the estimated T |
| `submit_experiments_5_cifar10.sh` | CIFAR-10 / CIFAR-10H application |
| `submit_experiments_6_bigearthnet.sh` | BigEarthNet application |

Each launcher selects one parameter block via a `CONF` variable (its leading comment names
the figure; see `experiments/metadata/map_full.md`). Override it on the command line:

```bash
cd experiments
CONF=3 ./submit_experiments_1.sh          # Figure 3 and Figures A19-A23
```

The launcher iterates over the parameter grid and submits one `sbatch` job per
configuration and random seed. Jobs already completed (an output `.txt` exists) are skipped.
Output `.txt` files are written under `experiments/results/expN/`.

**Running a single repetition locally** (no cluster): call the per-experiment wrapper
directly with the positional arguments the launcher would pass, e.g.

```bash
cd experiments
bash exp_1.sh 1 synthetic1 20 4 1.0 RFC 0.1 0.2 uniform 10000 1000 none 1
#             CONF DATA  p  K sig model eps  nu  contam   n_tr n_cal est seed
```

The full grid (all seeds and configurations, thousands of jobs) is not practical on a single
machine -- this is why Stage 1 targets a cluster.

**Data for the real-data experiments.** Before running `submit_experiments_5_cifar10.sh` or
`submit_experiments_6_bigearthnet.sh`, run the download scripts in `data/` (Step 0b).
BigEarthNet in particular is **> 30 GB** and should be fetched on a cluster.

## Step 1b -- Collect the results for plotting

`experiments/make_plots.R` reads from `experiments/results_hpc/`. On the cluster the jobs
write to `experiments/results/`; gather those into one `experiments/results_hpc/` tree.
For an `rsync` helper, copy `experiments/download.sh.example` to `experiments/download.sh`
and edit `REMOTE` (your `download.sh` is git-ignored).

## Step 2 -- Build the figures

Most figures are built by `experiments/make_plots.R`; the Richardson figures (A13-A14) are
built by `experiments/exp_richardson.ipynb` itself. See `experiments/metadata/map_full.md`
for the figure -> code mapping.

```bash
cd experiments
Rscript make_plots.R
```

runs top to bottom and writes every manuscript figure to `experiments/figures/` (PDF/PNG).
It reads `experiments/results_hpc/expN/` for most figures and
`experiments/results/simplified_methods/` for the finite-sample-correction figures A1-A6.
To rebuild a single figure, run just its section (use `map_full.md` to find it).

## Expected run time

- **Step 2 only** (rebuild all paper figures from a populated `results_hpc/`): a few minutes
  on a laptop.
- **Step 1 from scratch** (all experiments): not feasible on a desktop -- thousands of
  independent jobs across datasets, parameter grids and seeds; intended for a Slurm cluster.
