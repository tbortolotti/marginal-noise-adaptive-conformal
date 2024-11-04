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
        "0" : 1,
        "1": 1,
        "2": 1,
        "3": 1,
        "4": 1,
        "5": 1,
        "6": 1,
        "7": 1,
        "8": 1,
        "9": 1,
        "10": 1,
        "11": 0,
        "12": 0,
        "13": 0,
        "14": 1,
        "15": 1,
        "16": 1,
        "17": 1,
        "18": 1,
        "19": 1,
        "20": 1,
        "21": 1,
        "22": 1,
        "23": 1,
        "24": 1,
        "25": 1,
        "26": 1,
        "27": 1,
        "28": 1,
        "29": 1,
        "30": 1,
        "31": 1,
        "32": 1,
        "33": 1,
        "34": 1,
        "35": 1,
        "36": 1,
        "37": 1,
        "38": 1,
        "39": 1,
        "40": 1,
        "41": 1,
        "42": 1
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
