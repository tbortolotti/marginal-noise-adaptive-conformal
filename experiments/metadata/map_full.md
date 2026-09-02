<!-- Reproducibility map: manuscript display item -> code that produces it. -->
<!-- Appendix numbering follows the revised manuscript: A1-A6 finite-sample correction, -->
<!-- A7-A12 comparison figures, A13-A14 Richardson, A16-A31 the remaining experiments. -->
<!-- Figure A15 is a schematic (NN architecture) with no numerical content -- see end note. -->
<!-- If the journal numbering changes again, this file is the single place to update. -->

Experiments are launched from the `submit_experiments_*.sh` scripts in `experiments/` (each
selects one parameter block via a `CONF` variable whose block comment names the figure(s)),
or from the two Jupyter notebooks. `experiments/submit_paper_experiments.sh` runs every
`(launcher, CONF)` pair below in one go (`--list` / `--dry-run` supported). Raw per-repetition
output is written to `experiments/results/` and collected into `experiments/results_hpc/`
for plotting (`experiments/download.sh`). Figures are drawn by `experiments/make_plots.R`,
except the Richardson figures (A13-A14), which the notebook that produces their results also
plots.

The `results_hpc/` column below is where `make_plots.R` reads each figure's data. These
folders are **not** tracked in the repository (size) -- regenerate them by running the
experiments (Stage 1). The pre-computed results also accompany the paper's supplementary
material.

Launcher -> experiment driver:

| Launcher | Driver |
|---|---|
| `submit_experiments_1.sh` | `exp_1.py` |
| `submit_experiments_2_comparisons.sh` | `exp_2_comparisons.py` |
| `submit_experiments_3_T_estimation.sh` | `exp_T_estimation.py` |
| `submit_experiments_4.sh` | `exp_with_estimated_T.py` |
| `submit_experiments_5_cifar10.sh` | `exp_cifar.py` |
| `submit_experiments_6_bigearthnet.sh` | `exp_bigearthnet.py` |

## Main-text figures

| Fig. | Launcher | `CONF` | Results dir | `make_plots.R` function |
|:------|:-----------------------------------------|:-------|:-------------------|:------------------------|
| 1 | `submit_experiments_1.sh` | 1 | `results_hpc/exp1/` | `make_figure_1` |
| 2 | `submit_experiments_1.sh` | 2 | `results_hpc/exp2/` | `make_figure_2` |
| 3 | `submit_experiments_1.sh` | 3 | `results_hpc/exp3/` | `make_figure_3` |
| 4 | `submit_experiments_1.sh` | 4 | `results_hpc/exp4/` | `make_figure_4` |
| 5 | `submit_experiments_4.sh` | 401 | `results_hpc/exp401/` | `make_figure_401` |
| 6 | `submit_experiments_6_bigearthnet.sh` | 601 | `results_hpc/exp601/` | `make_figure_601b` |

## Appendix figures A1-A6 -- finite-sample correction term

Results produced by the notebook `experiments/exp_simplified_methods.ipynb` (writes
`experiments/results/simplified_methods/`). Figures drawn by the "Simplified Methods for
Special Contamination Models" section at the end of `experiments/make_plots.R`.

| Fig. | `make_plots.R` function | Output |
|:---|:---|:---|
| A1 | `make_figure_A1` | `figures/delta_marg_const_B_Klab_eps0.100000.pdf` |
| A2 | `make_figure_A2` | `figures/delta_marg_const_B_nlab_eps0.100000.pdf` |
| A3 | `make_figure_A3` | `figures/delta_marg_const_BRR_Klab_eps0.100000.pdf` |
| A4 | `make_figure_A4` | `figures/delta_marg_const_BRR_nlab_eps0.100000.pdf` |
| A5 | `make_figure_A5` (eps=0.1) | `figures/delta_marg_const_BRR_nu_var_eps0.100000.pdf` |
| A6 | `make_figure_A5` (eps=0.5) | `figures/delta_marg_const_BRR_nu_var_eps0.500000.pdf` |

## Appendix figures A7-A12 -- comparison with label-conditional / Clarkson et al. (2024)

| Fig. | Launcher | `CONF` | Results dir | `make_plots.R` function |
|:------|:-----------------------------------------|:-------|:-------------------|:------------------------|
| A7  | `submit_experiments_2_comparisons.sh` | 203 | `results_hpc/exp203/` | `make_figure_203` |
| A8  | `submit_experiments_2_comparisons.sh` | 204 | `results_hpc/exp204/` | `make_figure_204` |
| A9  | `submit_experiments_2_comparisons.sh` | 203 | `results_hpc/exp203/` | `make_figure_203` |
| A10 | `submit_experiments_2_comparisons.sh` | 204 | `results_hpc/exp204/` | `make_figure_204` |
| A11 | `submit_experiments_2_comparisons.sh` | 203 | `results_hpc/exp203/` | `make_figure_203` |
| A12 | `submit_experiments_2_comparisons.sh` | 204 | `results_hpc/exp204/` | `make_figure_204` |

## Appendix figures A13-A14 -- Richardson extrapolation

Both the results and the figures are produced by the notebook
`experiments/exp_richardson.ipynb` (cells "Figure A13" and "Figure A14"; writes
`experiments/results/richardson_extrapolation/`, saves the PNGs under `figures/`).

## Appendix figures A16-A31 -- remaining experiments

These figures keep their numbers from the earlier version of the manuscript.

| Fig. | Launcher | `CONF` | Results dir | `make_plots.R` function |
|:------|:-----------------------------------------|:-------|:-------------------|:------------------------|
| A16, A17 | `submit_experiments_1.sh` | 1 | `results_hpc/exp1/` | `make_figure_1` (uniform / block contamination) |
| A18 | `submit_experiments_1.sh` | 2 | `results_hpc/exp2/` | `make_figure_2A` |
| A19-A23 | `submit_experiments_1.sh` | 3 | `results_hpc/exp3/` | `make_figure_3` |
| A24 | `submit_experiments_2_comparisons.sh` | 201 | `results_hpc/exp201/` | `make_figure_201` |
| A25, A26 | `submit_experiments_4.sh` | 403 | `results_hpc/exp403/` | `make_figure_403` |
| A27 | `submit_experiments_4.sh` | 406 | `results_hpc/exp406/` | `make_figure_406` |
| A28 | `submit_experiments_3_T_estimation.sh` | 301 | `results_hpc/exp301/` | `make_figure_301` |
| A29 | `submit_experiments_3_T_estimation.sh` | 304 | `results_hpc/exp304/` | `make_figure_304` |
| A30 | `submit_experiments_3_T_estimation.sh` | 305 | `results_hpc/exp305/` | `make_figure_305` |
| A31 | `submit_experiments_5_cifar10.sh` | 503 | `results_hpc/exp503/` | `make_figure_503` |

## In-text numbers

The numerical values quoted in the text (coverage / prediction-set-size summaries) are read
off the same `results_hpc/expN/` folders and the same `make_plots.R` summaries used for the
corresponding figures; no separate computation is involved.

## Method implementation

The methods themselves are implemented in the `cln/` Python package:

| Method (paper name) | Implementation |
|---|---|
| Adaptive / Adaptive+ (optimized) | `cln/classification.py`, `cln/optimization.py` |
| Adaptive (simplified) | `cln/classification.py` |
| Adaptive (asymptotic) | `cln/asymptotic.py` |
| Label-conditional noise-adaptive | `cln/classification_label_conditional.py` |
| Contamination models (uniform / block / RRB / ...) | `cln/contamination.py` |
| Noise-transition-matrix estimation | `cln/T_estimation.py`, `cln/T_estimation_NN.py`, `cln/estimation.py` |
| Synthetic data generators (`synthetic1`-`synthetic6`) | `cln/data.py` |
| Standard split conformal + black boxes (baseline) | `third_party/arc/` |
| Clarkson et al. (2024) competitor | `third_party/clarkson/` |

## Tables

The manuscript contains **no code-generated tables** (confirmed by the authors).

## Not included in the manuscript

`CONF` blocks `0`, `200`, `300`, `400`, `500`, `600` are one-line smoke tests.
`CONF` blocks `202`, `302`, `303`, `402`, `404`, `405`, `501`, `502` are exploratory runs
that do not appear in the manuscript; their `make_figure_*` calls in `make_plots.R` are
commented out, so `Rscript make_plots.R` produces only the manuscript figures.

## Note -- Figure A15

Figure A15 is a schematic diagram of the neural-network architecture; it contains no
numerical results and is not produced by code.
