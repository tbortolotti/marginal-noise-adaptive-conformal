import logging
import hydra
import torch
import pytorch_lightning as pl
from hydra.utils import instantiate
from hydra.core.hydra_config import HydraConfig  # Import HydraConfig
from omegaconf import DictConfig
from pathlib import Path

from bigearthnet.datamodules.bigearthnet_datamodule import BigEarthNetDataModule
from bigearthnet.models.bigearthnet_module import BigEarthNetModule

logger = logging.getLogger(__name__)

@hydra.main(config_path="configs", config_name="config", version_base="1.2")
def main(cfg: DictConfig):
    logger.info("Evaluating model...")

    label_mapping = {
        "0" : "Agriculture",
        "1": "Other",
        "2": "Agriculture",
        "3": "Other",
        "4": "Coast, water and wetlands",
        "5": "Vegetation",
        "6": "Other",
        "7": "Coast, water and wetlands",
        "8": "Agriculture",
        "9": "Vegetation",
        "10": "Other",
        "11": "Urban fabric",
        "12": "Urban fabric",
        "13": "Other",
        "14": "Coast, water and wetlands",
        "15": "Agriculture",
        "16": "Other",
        "17": "Industrial or commercial units",
        "18": "Coast, water and wetlands",
        "19": "Other",
        "20": "Agriculture",
        "21": "Other",
        "22": "Vegetation",
        "23": "Vegetation",
        "24": "Vegetation",
        "25": "Arable land",
        "26": "Agriculture",
        "27": "Agriculture",
        "28": "Coast, water and wetlands",
        "29": "Arable land",
        "30": "Other",
        "31": "Arable land",
        "32": "Other",
        "33": "Coast, water and wetlands",
        "34": "Coast, water and wetlands",
        "35": "Vegetation",
        "36": "Coast, water and wetlands",
        "37": "Vegetation",
        "38": "Other",
        "39": "Vegetation",
        "40": "Agriculture",
        "41": "Coast, water and wetlands",
        "42": "Coast, water and wetlands"
        }

    # Load the model from the checkpoint
    model = BigEarthNetModule(cfg)
    model.load_state_dict(torch.load("trained_model.pth"))
    

    # fetch the transforms used in the model
    transforms = instantiate(cfg.transforms.obj)

    # instantiate the datamodule
    datamodule = BigEarthNetDataModule(
        cfg.datamodule.dataset_dir,
        cfg.datamodule.dataset_name,
        cfg.datamodule.batch_size,
        cfg.datamodule.num_workers,
        transforms,
        label_mapping
    )
    datamodule.setup()

    trainer = pl.Trainer(accelerator=cfg.trainer.accelerator, devices=cfg.trainer.devices)

    # Evaluate best model on test set
    trainer.test(model=model, datamodule=datamodule)
    logger.info("Test evaluation Done.")

if __name__ == "__main__":
    main()
