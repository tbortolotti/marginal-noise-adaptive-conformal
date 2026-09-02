<!-- Data sources and access. Referenced by README.md. -->

The paper uses three kinds of data. None of it is redistributed in this repository; all of
it is openly available and is obtained as described below. Random seeds are set in every
experiment script (integer `seed` argument, threaded through as `random_state`), so the
synthetic data and the label-noise draws are exactly reproducible.

> **Data availability statement (from the manuscript).** The BigEarthNet and CIFAR-10H data
> that support the findings of this study are openly available at <https://bigearth.net/> and
> <https://github.com/jcpeterson/cifar-10h>; see Sumbul et al. (2019) and Peterson et al. (2019).

Everything below the raw-image download is data *pre-processing* (label refinement, grouping,
classifier training) and is outside the scope of reproducing the paper's results: the
processed label files and the trained classifier are included in this repository.

## 0. Downloading the real data -- `data/`

The `data/` folder holds download scripts (only the scripts are tracked). From the repo root:

```bash
cd data
./download_all.sh                     # CIFAR-10 (~170 MB) + CIFAR-10H (~1 MB)
./download_all.sh --with-bigearthnet  # also BigEarthNet  (> 30 GB -- see warning below)
```

| Script | Into | Size |
|---|---|---|
| `data/download_cifar10.sh` | `data/cifar10/` | ~170 MB |
| `data/download_cifar10h.sh` | `data/cifar10h/cifar10h-counts.npy` | ~1 MB |
| `data/download_bigearthnet.sh --yes` | `data/bigearthnet-full/` | **> 30 GB** |

The experiment scripts read from `data/` by default (see the tables in sections 2 and 3).
None of this is needed for Stage 2 (building figures) once the experiment results are in
place. See `data/README.md` for requirements.

## 1. Synthetic data (`synthetic1`-`synthetic6`)

Fully generated in code -- no external files. The generative models are in `cln/data.py`
(`DataModel_1` ... `DataModel_6`); the noise-transition matrices and the contamination
process are in `cln/contamination.py`. Re-running with the same `seed` reproduces the same
dataset.

## 2. CIFAR-10 / CIFAR-10H (Appendix CIFAR experiment)

- **Images and clean labels -- CIFAR-10.** Standard dataset, fetched by
  `data/download_cifar10.sh` (torchvision) into `data/cifar10/cifar-10-batches-py/`.
- **Human annotation counts -- CIFAR-10H.** The file `cifar10h-counts.npy` (per-image
  human-label histograms) from <https://github.com/jcpeterson/cifar-10h>
  (Peterson, Battleday, Griffiths & Russakovsky, 2019), fetched by
  `data/download_cifar10h.sh` into `data/cifar10h/`. The BEN-1.0-style "noisy" labels are
  drawn from these histograms; the CIFAR-10 labels are the clean target.
- **Classifier.** Pretrained CIFAR-10 ResNet-18 from `cifar10_models`, vendored under
  `third_party/pytorch-cifar10/` (huyvnphan/PyTorch_CIFAR10, MIT; see its `PROVENANCE.md`).
  The package code and the `resnet18.pt` weights must be present at
  `third_party/pytorch-cifar10/cifar10_models/` -- if that folder is missing from your copy
  of this repository, add it as described in `PROVENANCE.md`. Override the location with the
  `PYTORCH_CIFAR10_PATH` environment variable.

**Where the data goes.** `experiments/exp_cifar.py` reads these two directories, defaulting
to the `data/` folder; override with environment variables:

| Variable | Default | Contents |
|---|---|---|
| `CIFAR10_DIR` | `data/cifar10` | `cifar-10-batches-py/` |
| `CIFAR10H_DIR` | `data/cifar10h` | `cifar10h-counts.npy` |

## 3. BigEarthNet, 6-class version (main real-data experiment, Figure 6)

### What the experiment consumes

The BigEarthNet experiment draws its calibration set -- and the sample used to estimate the
noise-transition matrix `T` -- from a pool of **~270,000 Sentinel-2 image patches**
(BigEarthNet-v1.0; Sumbul et al., 2019). Each 120x120-pixel patch covers a 1.2x1.2 km area
and carries two label sets:

- the **original BigEarthNet 1.0** CORINE labels, treated as the *noisy* labels, and
- the **Refined BigEarthNet (REBEN)** labels (Clasen et al., 2024), treated as the *clean*
  target.

Both label sets are reduced to a **6-class single-label** scheme
(*Coast, water and wetlands*; *Arable land*; *Agriculture*; *Vegetation*; *Urban fabric*;
plus *Mixed* for patches spanning categories), with severe class imbalance (frequencies from
~0.69 to ~0.001).

### How to obtain the data

1. **Patch lists (included).** The patch IDs of the ~270,000-patch pool are listed in
   `third_party/bigearthnet/data/scripts/splits/train.csv` (269,695 rows; the original
   BigEarthNet train split from TU Berlin). The `val.csv` / `test.csv` splits sit next to it.
2. **Per-patch labels (included).** `third_party/bigearthnet/data/train_bigearthnet-full.csv`
   (and `val_`/`test_`) give, for each patch, the noisy and clean grouped labels in columns
   `v1-labels-grouped` and `v2-labels-grouped` (0-5 = the six classes).
3. **19 -> 6 label mapping (included).**
   `third_party/bigearthnet/data/label_mapping.json` maps each original CORINE label (keys
   index into `third_party/bigearthnet/data/class_list.json`) to one of the six classes, or to
   `"Other"`. This is the file that defines the grouping.
4. **The image patches themselves (must be downloaded -- LARGER THAN 30 GB, use a cluster).**
   The experiment reads a pre-processed ("hub"-format) copy of the patches. Fetch it with:

   ```bash
   cd data
   ./download_bigearthnet.sh --yes       # gdown the archive from Google Drive
   ```

   (or, as a Slurm batch job, `sbatch third_party/bigearthnet/download_bigearthnet_full.sh`).
   This downloads and extracts the `bigearthnet-full` archive into `data/bigearthnet-full/`.
   The compressed archive alone is > 30 GB and the extracted data is larger still -- run it on
   a machine/cluster with ample free disk and a stable connection, **not a laptop**.
   `experiments/exp_bigearthnet.py` finds it via
   `third_party/bigearthnet/configs/config.yaml` (`datamodule.dataset_dir: ../data/`,
   `datamodule.dataset_name: bigearthnet-full`). `bigearthnet-mini` / `-medium` are much
   smaller subsets (same `GDRIVE_URLS` table) useful for a smoke test.
5. **Classifier (included).** The non-conformity scores use a ResNet-34 trained on 25,000
   patches with noisy labels (adapting the code of Pinto & St-Charles, 2022). Training is
   pre-processing and out of scope; the trained checkpoint
   `third_party/bigearthnet/trained_model.pth` is included and loaded directly by
   `experiments/exp_bigearthnet.py`.

The downloaded patches (`data/bigearthnet-full/`) are **not** stored in this repository.

## Data dictionary

Synthetic variables: documented in `cln/data.py`. CIFAR-10 / CIFAR-10H / BigEarthNet: see the
upstream documentation at the links above, the manuscript appendix, and the label files noted
in section 3.

## References for the data

- J. C. Peterson, R. M. Battleday, T. L. Griffiths, O. Russakovsky (2019). *Human uncertainty
  makes classification more robust.* Proc. IEEE/CVF ICCV, 9617-9626.
- G. Sumbul, M. Charfuelan, B. Demir, V. Markl (2019). *BigEarthNet: A large-scale benchmark
  archive for remote sensing image understanding.* IGARSS 2019, 5901-5904.
- A. Clasen et al. (2024). *Refined BigEarthNet (REBEN).*
- J. Pinto, P.-L. St-Charles (2022). *BigEarthNet.* <https://github.com/jerpint/bigearthnet>
