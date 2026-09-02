<!-- Software requirements. Referenced by README.md. -->

### Version of primary software used

- **Python**: 3.11 for the synthetic-data and CIFAR-10 experiments, the `cln/` package,
  and all method implementations. 3.8 for the BigEarthNet experiment (`third_party/bigearthnet`
  requires an older stack). Two conda environments are used, referenced by the
  `experiments/exp_*.sh` wrappers as `default` and `bigearth`.
- **R**: 4.4.2 for `experiments/make_plots.R` (all figures).

### Python -- synthetic + CIFAR-10 (`default` environment)

Used by `cln/`, `experiments/exp_1.py`, `exp_2_comparisons.py`, `exp_T_estimation.py`,
`exp_with_estimated_T.py`, `exp_cifar.py`, and the two notebooks.

| Package | Version | Paper minimum (from the manuscript) |
|---|---|---|
| numpy | 1.26.4 | >= 1.25.0 |
| scipy | 1.14.0 | >= 1.11.1 |
| scikit-learn | 1.5.1 | >= 1.3.0 |
| pandas | 2.2.2 | >= 2.0.3 |
| statsmodels | 0.14.2 | >= 0.14.0 |
| cvxpy | 1.5.2 | >= 1.5.2 |
| tqdm | 4.66.5 | >= 4.65.0 |
| matplotlib | 3.9.1 | -- |
| torch | 2.4.0 | >= 1.10.2 |
| torchvision | 0.19.0 | (matched to torch) |

`torch` / `torchvision` are only needed for the CIFAR-10 experiment and the neural
noise-matrix estimator (`cln/T_estimation_NN.py`); the synthetic-data experiments run
without them. The CIFAR-10 classifier is the pretrained ResNet-18 from `cifar10_models`,
vendored under `third_party/pytorch-cifar10/` (huyvnphan/PyTorch_CIFAR10, MIT; see
`third_party/pytorch-cifar10/PROVENANCE.md` -- the package code + `resnet18.pt` weights must
be present there). Loading its state dict works with recent `torch`; upstream targets
`torch>=1.7`.

Install with `pip install -r requirements.txt` (repo root). `seaborn` is used only by
`experiments/exp_richardson.ipynb`.

### Python -- BigEarthNet (`bigearth` environment)

Used by `experiments/exp_bigearthnet.py` and `third_party/bigearthnet/`.

Versions are those of the conda environment that produced the paper's BigEarthNet results.

| Package | Version |
|---|---|
| python | 3.8 |
| numpy | 1.24.4 |
| scipy | 1.10.1 |
| scikit-learn | 1.3.2 |
| pandas | 2.0.3 |
| torch | 2.0.1 |
| torchvision | 0.15.2 |
| pytorch-lightning | 2.0.0 |
| timm | 0.6.12 |
| torchgeo | 0.4.1 (imported by `third_party/bigearthnet/models/bigearthnet_module.py`) |
| kornia | 0.6.12 (imported by the same module) |
| hydra-core | 1.3.2 |
| omegaconf | 2.3.0 |
| gdown | 5.2.0 |
| hub | 3.0.1 (activeloop "hub"; reads the pre-processed BigEarthNet dataset) |

Install with `pip install -r requirements-bigearthnet.txt` (repo root).

### R (figures)

| Package | Version |
|---|---|
| tidyverse | 2.0.0 (ggplot2 4.0.1, dplyr 1.1.4, tidyr 1.3.2, readr 2.1.5) |
| latex2exp | 0.9.8 |
| RColorBrewer | 1.1.3 |
| cowplot | 1.1.3 |
| patchwork | 1.3.0 |
| ggforce | 0.5.0 |

A full list is in `metadata/R-installed-packages.csv`; `sessionInfo()` output is in
`metadata/R-sessionInfo.txt`.

### Interfacing R and Python

The synthetic-data and CIFAR experiments are pure Python; figures are pure R and read
plain-text (`.txt`, comma-delimited) result files. There is **no** in-process R/Python
bridge (no `rpy2`) in the paper workflow.

### External / non-CRAN components (vendored under `third_party/`)

- `third_party/arc/` -- split-conformal baselines and black-box wrappers.
- `third_party/clarkson/` -- implementation of the Clarkson et al. (2024) competitor.
- `third_party/pytorch-cifar10/` -- huyvnphan/PyTorch_CIFAR10 (MIT), the pretrained CIFAR-10
  ResNet-18 used by the CIFAR appendix experiment. See its `PROVENANCE.md`; the
  `cifar10_models/` package + `resnet18.pt` weights must be added there.
- `third_party/bigearthnet/` -- BigEarthNet data module and model, modified for the
  6-class setting, with the trained checkpoint `trained_model.pth`.
