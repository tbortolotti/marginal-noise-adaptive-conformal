# `data/` -- real datasets

The scripts here download the real datasets into this folder. Nothing here is
tracked in the repository except the scripts themselves. See `../metadata/data.md`
for full sourcing and citations.

| Script | Fetches | Into | Size |
|---|---|---|---|
| `download_cifar10.sh` | CIFAR-10 (torchvision) | `data/cifar10/` | ~170 MB |
| `download_cifar10h.sh` | CIFAR-10H human-annotation counts | `data/cifar10h/cifar10h-counts.npy` | ~1 MB |
| `download_bigearthnet.sh` | pre-processed BigEarthNet ("hub" format) | `data/bigearthnet-full/` | **> 30 GB** |
| `download_all.sh` | CIFAR-10 + CIFAR-10H (add `--with-bigearthnet` for BEN) | as above | |

## Usage

```bash
cd data
./download_all.sh                       # CIFAR-10 + CIFAR-10H only
./download_all.sh --with-bigearthnet     # also BigEarthNet (> 30 GB)
# or individually:
./download_bigearthnet.sh --yes
```

## Requirements

- `download_cifar10.sh` needs a Python with `torch` + `torchvision` (the `default`
  environment).
- `download_cifar10h.sh` needs `curl`.
- `download_bigearthnet.sh` needs `gdown` (`pip install gdown`; already in the
  `bigearth` environment).

## Notes

- **BigEarthNet is larger than 30 GB.** Run its download on a cluster or a machine
  with ample free disk and a stable connection, not a laptop. It is only needed to
  re-run the BigEarthNet experiment (`experiments/exp_bigearthnet.py`); the paper
  figure is rebuilt from the shipped results without it.
- The experiment scripts read from here by default:
  `experiments/exp_cifar.py` uses `data/cifar10` and `data/cifar10h` (override with
  `CIFAR10_DIR` / `CIFAR10H_DIR`); `third_party/bigearthnet/configs/config.yaml`
  points `datamodule.dataset_dir` at `../data`.
- The BigEarthNet Google Drive id used by `download_bigearthnet.sh` is the
  `bigearthnet-full` entry of `GDRIVE_URLS` in
  `../third_party/bigearthnet/datamodules/bigearthnet_datamodule.py`.
