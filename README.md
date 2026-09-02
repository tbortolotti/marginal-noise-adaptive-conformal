# Noise-Adaptive Conformal Classification with Marginal Coverage

Software implementation of the methods and experiments from the paper:

> **Noise-Adaptive Conformal Classification with Marginal Coverage**
> T. Bortolotti, Y. X. Wang, X. Tong, A. Menafoglio, S. Vantini, M. Sesia
> arXiv: <https://arxiv.org/abs/2501.18060>

## Abstract

Conformal inference provides a rigorous statistical framework for uncertainty quantification
in machine learning, enabling well-calibrated prediction sets with precise coverage
guarantees for any classification model. However, its reliance on the idealized assumption of
perfect data exchangeability limits its effectiveness in the presence of real-world
complications, such as low-quality labels -- a widespread issue in modern large-scale data
sets. This work tackles this open problem by introducing an adaptive conformal inference
method capable of efficiently handling deviations from exchangeability caused by random label
noise, leading to informative prediction sets with tight marginal coverage guarantees even in
those challenging scenarios. We validate our method through extensive numerical experiments
demonstrating its effectiveness on synthetic and real data sets, including CIFAR-10H and
BigEarthNet.

## Repository layout

| Path | Contents |
|---|---|
| `cln/` | Python package implementing the noise-adaptive conformal methods |
| `data/` | Scripts to download the real datasets (CIFAR-10, CIFAR-10H, BigEarthNet) |
| `experiments/` | Experiment drivers, Slurm launchers, and the plotting script `make_plots.R` |
| `experiments/metadata/map_full.md` | Reproducibility map: each figure -> the code that makes it |
| `metadata/` | Dependency lists, workflow description, data sourcing |
| `requirements.txt`, `requirements-bigearthnet.txt` | Python dependencies for the two environments |
| `third_party/` | Vendored dependencies (`arc`, `clarkson`, `pytorch-cifar10`, `bigearthnet`) |
| `LICENSE` | MIT License |

## Prerequisites

See [`metadata/dependencies.md`](metadata/dependencies.md). Summary: Python 3.11
(`pip install -r requirements.txt`) plus a separate Python 3.8 environment for the
BigEarthNet experiment (`requirements-bigearthnet.txt`), and R 4.4. R package snapshots are
in `metadata/`.

## Data

See [`metadata/data.md`](metadata/data.md). The synthetic data is generated in code (seeded).
CIFAR-10 / CIFAR-10H and BigEarthNet are publicly available; none of it is redistributed
here. Download scripts are in [`data/`](data/README.md):

```bash
cd data
./download_all.sh                     # CIFAR-10 + CIFAR-10H
./download_all.sh --with-bigearthnet  # BigEarthNet -- larger than 30 GB, use a cluster
```

The BigEarthNet classifier checkpoint (`third_party/bigearthnet/trained_model.pth`) *is*
included, so that experiment does not need the classifier retrained.

## Reproducibility workflow

Full details in [`metadata/workflow.md`](metadata/workflow.md). In brief:

- **Step 1 -- run experiments** (Slurm cluster). Each family of experiments has a
  `experiments/submit_experiments_*.sh` launcher with a `CONF` switch that selects the
  configuration for a given figure; `experiments/submit_paper_experiments.sh` runs them all.
  Outputs go to `experiments/results/`; collect them for plotting into
  `experiments/results_hpc/` (`experiments/download.sh`). Not feasible on a single desktop.
- **Step 2 -- build figures** (laptop, minutes once `results_hpc/` is populated). Run
  `experiments/make_plots.R`; it reads `experiments/results_hpc/` (the finite-sample-correction
  figures A1--A6 read `experiments/results/simplified_methods/`). Use the reproducibility map
  to find the section for a given figure.

> The experiment outputs (`experiments/results/`, `experiments/results_hpc/`) are **not**
> tracked in this repository, because of their size. Regenerate them with Step 1. The
> pre-computed results also accompany the paper's supplementary material.

## Reproducibility map

The mapping from each manuscript figure to the code that produces it is in
[`experiments/metadata/map_full.md`](experiments/metadata/map_full.md).

## Expected run time

- Step 2 only, once `results_hpc/` is populated: a few minutes.
- Step 1 from scratch (all experiments): not feasible on a desktop; intended for a cluster.

## License

MIT License -- see [`LICENSE`](LICENSE).
