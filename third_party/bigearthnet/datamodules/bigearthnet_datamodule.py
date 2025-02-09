import logging
import os
import pathlib
import tarfile
import typing

import gdown
import hub
import numpy as np
import random
import torch
import pytorch_lightning as pl
import torch.utils.data.dataloader
import torch.utils.data.dataset
from torch.utils.data import SequentialSampler, SubsetRandomSampler

logger = logging.getLogger(__name__)

GDRIVE_URLS = {
    "bigearthnet-mini": "https://drive.google.com/file/d/1X2-NpZ4ExUooAi8tkCBWAlmHZAG_N3Ux/view?usp=sharing",
    "bigearthnet-medium": "https://drive.google.com/file/d/1YW4ugRQTl-YF_ZpLO7gIRSlHnB2Cwslz/view?usp=sharing",
    "bigearthnet-full": "https://drive.google.com/file/d/1isUcPQvCn1xc5GWEmPDOqj_l24QH2ZtA/view?usp=sharing",
}

def download_data(dataset_dir, dataset_name):
    """Download and extract the specified dataset to dataset_dir if not already present."""
    assert dataset_name in GDRIVE_URLS.keys(), "dataset_name isn't available."
    dataset_path = pathlib.Path(os.path.join(dataset_dir, dataset_name))
    tar_fname = str(dataset_path.resolve()) + ".tar"

    # Skip download if dataset was already extracted
    if os.path.isdir(dataset_path):
        logger.info(f"Dataset already present at {dataset_path}, skipping download.")
        return dataset_path

    # Create dataset dir if it doesn't already exist
    if not os.path.isdir(dataset_dir):
        os.mkdir(dataset_dir)

    # download the dataset from google drive
    if not os.path.isfile(tar_fname):
        logger.info(f"Downloading {dataset_name} dataset to {tar_fname}")
        url = GDRIVE_URLS[dataset_name]
        gdown.download(url, str(tar_fname), fuzzy=True)

    # extract tar file
    with tarfile.open(str(tar_fname), "r") as tar:
        tar.extractall(path=str(dataset_dir), members=tar)

    logger.info(
        f"Succesfully downloaded and extracted {dataset_name} to {dataset_path}"
    )

    return dataset_path


def map_labels_to_single(hub_labels, mapping_dict):
    """Map multi-labels to a single label based on a given mapping.
    If multi_labels contains only labels that are mapped to 0, return 0.
    If multi_labels contains only labels that are mapped to 1, return 2.
    If multi_labels contains a mix of labels mapped to 0 and 1, return 1.
    """
    # Flatten hub_labels if it contains numpy arrays
    flat_labels = []
    for label in hub_labels:
        if isinstance(label, np.ndarray):
            flat_labels.extend(str(lbl) for lbl in label.tolist())
        else:   
            flat_labels.append(str(label))
    mapped_labels = {mapping_dict[label] for label in flat_labels if label in mapping_dict and mapping_dict[label]!="Other"}
    
    if mapped_labels=={"Coast, water and wetlands"}:
        return 0
    elif mapped_labels=={"Arable land"}:
        return 1
    elif mapped_labels=={"Agriculture"}:
        return 2
    elif mapped_labels=={"Vegetation"}:
        return 3
    elif len(mapped_labels)==1:
        return 4 # Other
    else:
        return 5 # Mix of two or more

class BigEarthNetHubDataset(torch.utils.data.dataset.Dataset):
    """Dataset class used to iterate over the BigEarthNet-S2 data."""

    def __init__(
        self,
        dataset_path: pathlib.Path,
        label_mapping: dict,
        transforms=None,
        **extra_hub_kwargs,
    ):
        """Initialize the BigEarthNet-S2 hub dataset (in read-only mode)."""
        self.dataset_path = dataset_path
        self.dataset = hub.load(self.dataset_path, read_only=True, **extra_hub_kwargs)
        self.transforms = transforms
        self.label_mapping = label_mapping

    def __len__(self) -> int:
        """Returns the total size (patch count) of the dataset."""
        return len(self.dataset)

    def __getitem__(self, idx: int) -> typing.Dict[str, typing.Any]:
        """Returns a single data sample loaded from the dataset.
        For BigEarthNet, the data sample simply consists of the patch data and its labels. The
        patch data and labels will be converted from their original types to float32 and int16,
        respectively, in order to make sure that PyTorch will be able to batch them.
        Labels are converted to one-hot representation.
        """
        item = self.dataset[
            int(idx)
        ]  # cast in case we're using numpy ints or something similar
        """
        # Print available attributes and methods for inspection
        print(dir(item))
    
        # Attempt to print specific attributes that may contain patch information
        if hasattr(item, 'info'):
            print("Info attribute:", item.info)
        if hasattr(item, 'patch_name'):
            print("Patch name:", item.patch_name)
        print(self.dataset.tensors)  # Check available tensors
        """
        assert tuple(self.tensor_names) == ("data", "labels")

        hub_labels = item["labels"].numpy()
        single_label = map_labels_to_single(hub_labels, self.label_mapping)
        labels = torch.tensor(single_label)

        img_data = item["data"].numpy().astype(np.float32)
        img_data = torch.tensor(img_data)

        if self.transforms:
            img_data = self.transforms(img_data)

        return {
            "data": img_data,
            "labels": labels,
        }

    def summary(self) -> None:
        """Forwards the call to print a summary of the dataset."""
        return self.dataset.summary()

    def visualize(self, *args, **kwargs):
        """Forwards the call to show the dataset content (notebook-only)."""
        return self.dataset.visualize(*args, **kwargs)

    @property
    def dataset_info(self) -> typing.Dict[str, typing.Any]:
        """Returns metadata information parsed from the hub dataset object."""
        return dict(self.dataset.info)

    @property
    def dataset_name(self) -> typing.AnyStr:
        """Returns the dataset name used to identify this particular dataset."""
        return self.dataset_info["name"]

    @property
    def class_names(self) -> typing.List[str]:
        """Returns the list of class names that correspond to the label indices in the dataset."""
        return list(self.dataset.info.class_names)

    @property
    def tensor_info(self) -> typing.Dict[str, typing.Dict[str, typing.Any]]:
        """Returns the dictionary of tensor info objects (hub-defined) parsed from the dataset.

        The returned objects can help downstream processing stages figure out what kind of data
        they will be receiving from this parser.
        """
        return {k: v.info for k, v in self.dataset.tensors.items()}

    @property
    def tensor_names(self) -> typing.List[str]:
        """Names of the tensors that will be provided in the loaded data samples.

        Note that additional tensors and other attributes may be loaded as well, but these are the
        'primary' fields that should be expected by downstream processing stages.
        """
        return list(self.tensor_info.keys())


class BigEarthNetDataModule(pl.LightningDataModule):
    """Data module class that prepares BigEarthNet-S2 dataset parsers and instantiates data loaders."""

    def __init__(
        self,
        dataset_dir: str,
        dataset_name: str,
        batch_size: int,
        num_workers: int = 0,
        transforms=None,
        label_mapping: dict = None,
        random_seed=None
    ):
        """Validates the hyperparameter config dictionary and sets up internal attributes."""
        super().__init__()
        self.dataset_name = dataset_name
        self.dataset_dir = pathlib.Path(dataset_dir)
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.train_dataset, self.valid_dataset, self.test_dataset = None, None, None
        self.transforms = transforms
        self.label_mapping = label_mapping
        self.random_seed=random_seed

    def setup(self, stage=None) -> None:
        """Parses and splits all samples across the train/valid/test datasets."""
        self.dataset_path = download_data(self.dataset_dir, self.dataset_name)

        if self.train_dataset is None:
            self.train_dataset = BigEarthNetHubDataset(
                self.dataset_path / "train",
                label_mapping=self.label_mapping,
                transforms=self.transforms,
            )
        if self.valid_dataset is None:
            self.valid_dataset = BigEarthNetHubDataset(
                self.dataset_path / "val",
                label_mapping=self.label_mapping,
                transforms=self.transforms,
            )
        if self.test_dataset is None:
            self.test_dataset = BigEarthNetHubDataset(
                self.dataset_path / "test",
                label_mapping=self.label_mapping,
                transforms=self.transforms,
            )

    def create_sampler(self, dataset_size: int, shuffle: bool = False):
        """Create a sampler for shuffling or deterministic sampling."""
        indices = list(range(dataset_size))
        if shuffle:
            generator = torch.Generator().manual_seed(self.random_seed)
            return SubsetRandomSampler(indices=indices, generator=generator)
        else:
            return SequentialSampler(indices)

    def train_dataloader(self) -> torch.utils.data.dataloader.DataLoader:
        """Creates the training dataloader using the training dataset."""
        assert self.train_dataset is not None, "must call 'setup' first!"

        # Create the sampler to shuffle the data
        sampler = self.create_sampler(len(self.train_dataset), shuffle=True)

        return torch.utils.data.dataloader.DataLoader(
            dataset=self.train_dataset,
            batch_size=self.batch_size,
            #batch_size=int(0.7*self.batch_size),
            #shuffle=True,
            sampler=sampler,
            num_workers=self.num_workers
        )

    def val_dataloader(self) -> torch.utils.data.dataloader.DataLoader:
        """Creates the validation dataloader using the validation data parser."""
        assert self.valid_dataset is not None, "must call 'setup' first!"

        # Create the sampler to shuffle the data
        sampler = self.create_sampler(len(self.train_dataset), shuffle=True)

        return torch.utils.data.dataloader.DataLoader(
            dataset=self.valid_dataset,
            #batch_size=self.batch_size,
            batch_size=int(0.15*self.batch_size),
            #shuffle=False,
            sampler=sampler,
            num_workers=self.num_workers,
        )

    def test_dataloader(self) -> torch.utils.data.dataloader.DataLoader:
        """Creates the testing dataloader using the testing data dataset."""
        assert self.test_dataset is not None, "must call 'setup' first!"

        # Create the sampler to shuffle the data
        sampler = self.create_sampler(len(self.train_dataset), shuffle=True)

        return torch.utils.data.dataloader.DataLoader(
            dataset=self.test_dataset,
            #batch_size=self.batch_size,
            batch_size=int(0.15*self.batch_size),
            #shuffle=False,
            sampler=sampler,
            num_workers=self.num_workers,
        )