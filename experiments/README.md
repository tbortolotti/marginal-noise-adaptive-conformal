# `experiments/` -- reproducing the manuscript figures

This directory contains everything needed to (1) run the experiments and (2) build the
figures. See `../metadata/workflow.md` for the full workflow and `metadata/map_full.md` for
the figure -> code map. Run all scripts from inside this directory.

## Running everything

`submit_paper_experiments.sh` submits **every experiment that backs a manuscript figure**
(all six launchers, the right `CONF` for each). It excludes the two notebook studies
(figures A1-A6, A13-A14). `--dry-run` prints without submitting; `--list` shows the
figure -> (launcher, `CONF`) table.

## Experiment drivers

| Driver | Wrapper | Slurm launcher | What it runs |
|---|---|---|---|
| `exp_1.py` | `exp_1.sh` | `submit_experiments_1.sh` | synthetic data: contamination strength / model / number of classes / marginal vs. label-conditional target |
| `exp_2_comparisons.py` | `exp_2_comparisons.sh` | `submit_experiments_2_comparisons.sh` | comparison with the label-conditional method and Clarkson et al. (2024) |
| `exp_T_estimation.py` | `exp_T_estimation.sh` | `submit_experiments_3_T_estimation.sh` | estimation of the noise-transition matrix T |
| `exp_with_estimated_T.py` | `exp_with_estimated_T.sh` | `submit_experiments_4.sh` | adaptive method using the estimated T |
| `exp_cifar.py` | `exp_cifar.sh` | `submit_experiments_5_cifar10.sh` | CIFAR-10 / CIFAR-10H application |
| `exp_bigearthnet.py` | `exp_bigearthnet.sh` | `submit_experiments_6_bigearthnet.sh` | BigEarthNet (6-class) application |

Each `submit_experiments_*.sh` selects one parameter block via a `CONF` variable whose
leading comment names the figure(s) it produces. Override it on the command line, e.g.
`CONF=3 ./submit_experiments_1.sh`. Each launcher submits one `sbatch` job per parameter
configuration and seed, skipping configurations whose output `.txt` already exists.

To run one repetition locally without Slurm, call the `.sh` wrapper directly with the
positional arguments the launcher passes (see `../metadata/workflow.md` for an example).

## Notebooks

- `exp_simplified_methods.ipynb` -- finite-sample-correction-term study (appendix figures
  **A1-A6**). Produces results only, into `results/simplified_methods/`; the figures
  are drawn by the "Simplified Methods" section at the end of `make_plots.R`.
- `exp_richardson.ipynb` -- Richardson-extrapolation study (appendix figures **A13-A14**).
  Produces both the results (`results/richardson_extrapolation/`) and the figures
  (PNGs under `figures/`).

## Helper scripts

- `data_torch.py` -- CIFAR-10 / CIFAR-10H dataset handling (used by `exp_cifar.py`); reads
  `../data/cifar10` and `../data/cifar10h` by default (see `../data/README.md`), overridable
  with `CIFAR10_DIR` / `CIFAR10H_DIR`.

## Plotting

- `make_plots.R` -- builds every manuscript figure. `Rscript make_plots.R` runs top to
  bottom; it reads `results_hpc/expN/` (and `results/simplified_methods/` for figures A1-A6)
  and writes to `figures/`. To rebuild one figure, run just its section -- use
  `metadata/map_full.md` to find it.

## Results and figures folders

None of these are tracked in the repository (size). Regenerate them by running the
experiments; the pre-computed results also accompany the paper's supplementary material.

- `results/` -- where the experiment jobs write:
  - `results/expN/` -- per-repetition outputs from the `submit_experiments_*.sh` launchers.
  - `results/simplified_methods/` -- outputs of `exp_simplified_methods.ipynb` (figs A1-A6).
  - `results/richardson_extrapolation/` -- outputs of `exp_richardson.ipynb` (figs A13-A14).
- `results_hpc/` -- `results/expN/` collected here for plotting (`download.sh`);
  `make_plots.R` reads `results_hpc/expN/` (and `results/simplified_methods/` for A1-A6).
- `figures/` -- output directory for `make_plots.R` (created on demand).
