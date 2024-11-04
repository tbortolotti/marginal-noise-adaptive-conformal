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
