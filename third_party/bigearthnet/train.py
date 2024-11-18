import logging

import hydra
import torch
import pytorch_lightning as pl
from hydra.utils import instantiate
from omegaconf import DictConfig

from bigearthnet.models.bigearthnet_module import BigEarthNetModule
from bigearthnet.datamodules.bigearthnet_datamodule import BigEarthNetDataModule
from bigearthnet.utils.reproducibility_utils import set_seed

log = logging.getLogger(__name__)


@hydra.main(config_path="configs", config_name="config", version_base="1.2")
def main(cfg: DictConfig):
    log.info("Beginning training...")

    # set seed if specified
    if cfg.experiment.get("seed"):
        set_seed(cfg.experiment.seed)

    # Define the label mapping dictionary
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

    # instantiate all objects from hydra configs
    model = BigEarthNetModule(cfg)

    # Modify the datamodule instantiation to pass the label mapping
    datamodule = BigEarthNetDataModule(
        dataset_dir=cfg.datamodule.dataset_dir,
        dataset_name=cfg.datamodule.dataset_name,
        batch_size=cfg.datamodule.batch_size,
        num_workers=cfg.datamodule.num_workers,
        transforms=instantiate(cfg.datamodule.transforms),
        label_mapping=label_mapping  # Pass the label mapping
    )
    
    datamodule.setup()

    # Manually create the Trainer object with the desired parameters
    trainer = pl.Trainer(
        max_epochs=cfg.trainer.max_epochs,
        accelerator=cfg.trainer.accelerator,
        devices=cfg.trainer.devices,
        num_sanity_val_steps=0,
        log_every_n_steps=cfg.trainer.log_every_n_steps,
    )

    # do the training
    trainer.fit(model, datamodule=datamodule)
    torch.save(model.state_dict(), "trained_model.pth")
    log.info("Training Done.")


if __name__ == "__main__":
    main()
