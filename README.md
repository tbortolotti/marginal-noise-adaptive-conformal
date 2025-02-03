# Marginal Adaptive Conformal Classification with Noisy Labels

This repository provides the software implementation of methods from the paper:

>  Noise-Adaptive Conformal Classification with Marginal Coverage
>  T. Bortolott, Y. X. Wang, X. Tong, A. Menafoglio, S. Vantini, M. Sesia
>  [arxiv preprint]{https://arxiv.org/abs/2501.18060}


## Abstract

Conformal inference provides a rigorous statistical framework for uncertainty quantification in machine learning, enabling well-calibrated prediction sets with precise coverage guarantees for any classification model. However, its reliance on the idealized assumption of perfect data exchangeability limits its effectiveness in the presence of real-world complications, such as low-quality labels -- a widespread issue in modern large-scale data sets. This work tackles this open problem by introducing an adaptive conformal inference method capable of efficiently handling deviations from exchangeability caused by random label noise, leading to informative prediction sets with tight marginal coverage guarantees even in those challenging scenarios. We validate our method through extensive numerical experiments demonstrating its effectiveness on synthetic and real data sets, including CIFAR-10H and BigEarthNet.

## Contents

- `cln/`: Python package implementing the marginal adaptive methods from the paper.
- `experiments/`: Code for reproducing the numerical experiments with simulated data.
- `third_party/`: Third-party Python packages used by this package, including the modified package for BigEarthNet analysis with six classes

## Prerequisites

For replicating the `experiments`:
- `numpy` (>= 1.25.0)
- `scipy` (>= 1.11.1)
- `scikit-learn` (>= 1.3.0)
- `pandas` (>= 2.0.3)
- `torch` (>= 1.10.2)
- `tqdm` (>= 4.65.0)
- `statsmodels` (>= 0.14.0)
- `cvxpy` (>= 1.5.2)

## Installation

Clone the development version from GitHub:

    git clone https://github.com/tbortolotti/marginal-noise-adaptive-conformal.git

## Reproducibility Instructions

See the [experiments/README.md](experiments/README.md) file for instructions on how to reproduce the paper's figures.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
