# Provenance -- `third_party/pytorch-cifar10/`

This folder is a vendored copy of **huyvnphan/PyTorch_CIFAR10**
(<https://github.com/huyvnphan/PyTorch_CIFAR10>), MIT License, (c) 2019 Huy Phan.
The upstream `LICENSE` and `README.md` are kept alongside this file.

## What it is used for

`experiments/data_torch.py` imports the pretrained CIFAR-10 ResNet-18 from this package:

```python
from cifar10_models.resnet import resnet18 as cifar_resnet18
```

It is used only for the **CIFAR-10 / CIFAR-10H appendix experiment**
(`experiments/exp_cifar.py`). It is not needed to rebuild any figure from the shipped
results.

## Contents

| Path | Origin |
|---|---|
| `README.md`, `LICENSE` | upstream, verbatim |
| `data.py`, `module.py`, `schduler.py`, `train.py` | upstream top-level files (not used by this project; kept for completeness) |
| `cifar10_models/` | upstream package: model definitions (`resnet.py`, ...) |
| `cifar10_models/state_dicts/resnet18.pt` | upstream pretrained weights for ResNet-18 (~43 MB) |

Upstream commit vendored: `641cac24371b17052b9bb6e56af1c83b5e97cd7f`.

Only the ResNet-18 weights are included. The other architectures' weights (VGG, DenseNet,
...) are not used and are omitted; obtain them from upstream
(`python train.py --download_weights 1`) if needed.

## Overriding the location

`data_torch.py` looks for this package at `../third_party/pytorch-cifar10` by default; set
the `PYTORCH_CIFAR10_PATH` environment variable to point elsewhere (e.g. to your own clone).
