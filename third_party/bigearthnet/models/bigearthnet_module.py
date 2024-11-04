import importlib.resources
import json
import logging
import os
import typing
import numpy as np
import copy

import pytorch_lightning as pl
import torch
from hydra.utils import instantiate
from omegaconf import DictConfig
from sklearn.metrics import (
    classification_report,
    multilabel_confusion_matrix,
    precision_recall_fscore_support,
)
from torch import optim


log = logging.getLogger(__name__)


class BigEarthNetModule(pl.LightningModule):
    """Base class for Pytorch Lightning model."""

    def __init__(self, cfg: DictConfig):
        super().__init__()
        self.cfg = cfg
        self.model = instantiate(cfg.model)

        # Identify the last fully connected layer
        if hasattr(self.model, 'fc'):  # For models with an 'fc' layer (e.g., ResNet)
            last_layer_size = self.model.fc.in_features
            self.model.fc = torch.nn.Linear(last_layer_size, 3)  # Adjust output to 3 classes
        elif hasattr(self.model, 'classifier'):  # For models with a 'classifier' layer
            last_layer_size = self.model.classifier.in_features
            self.model.classifier = torch.nn.Linear(last_layer_size, 3)  # Adjust output to 3 classes
        else:
            raise AttributeError("The model does not have an 'fc' or 'classifier' layer.")
        
        self.save_hyperparameters(cfg, logger=False)
        self.init_loss()

        self.validation_step_outputs = []
        self.test_step_outputs = []

    def init_loss(self):
        weights_file = self.cfg.loss.get("class_weights")
        if weights_file:
            # If specified in the config, the loss will be rebalanced according to data.
            assert os.path.isfile(weights_file)
            log.info(f"loading {weights_file} as weights for the loss function")
            with open(weights_file, "r") as f:
                data = json.load(f)
            pos_weight = torch.tensor(list(data.values()))
        else:
            pos_weight = None
        # self.loss_fn = torch.nn.BCEWithLogitsLoss(pos_weight=pos_weight)
        self.loss_fn = torch.nn.CrossEntropyLoss()

    def configure_optimizers(self):
        name = self.cfg.optimizer.name
        lr = self.cfg.optimizer.lr
        if name == "adam":
            optimizer = optim.Adam(
                self.model.parameters(),
                lr=lr,
            )
        elif name == "sgd":
            optimizer = optim.SGD(self.model.parameters(), lr=lr)
        else:
            raise ValueError(f"optimizer {name} not supported")
        return optimizer

    def _generic_step(self, batch, batch_idx):
        """Runs the prediction + evaluation step for training/validation/testing."""
        inputs = batch["data"]
        targets = batch["labels"]
        logits = self.model(inputs)
        #loss = self.loss_fn(logits, targets.float())
        loss = self.loss_fn(logits, targets)
        return {"loss": loss, "targets": targets, "logits": logits}

    def _generic_epoch_end(self, step_outputs):
        all_targets = []
        all_preds = []
        all_loss = []
        for outputs in step_outputs:
            logits = outputs["logits"]
            targets = outputs["targets"]
            #preds = torch.sigmoid(logits) > 0.5
            preds = torch.argmax(logits, dim=1)
            all_targets.extend(targets.cpu().numpy())
            #all_preds.extend(preds.type(targets.dtype).cpu().numpy())
            all_preds.extend(preds.cpu().numpy())

            loss = outputs["loss"]
            all_loss.append(loss.cpu().numpy())

        accuracy = (np.array(all_preds) == np.array(all_targets)).mean()
        avg_loss = sum(all_loss) / len(all_loss)

        metrics = {
            "accuracy": accuracy,
            "loss":avg_loss,
        }
        return metrics

    def training_step(self, batch, batch_idx):
        """Runs a prediction step for training, returning the loss."""
        outputs = self._generic_step(batch, batch_idx)
        self.log(
            "loss/train",
            outputs["loss"],
            on_step=True,
            on_epoch=True,
            prog_bar=True,
            logger=True,
        )
        return outputs

    def on_training_epoch_end(self, training_step_outputs):
        metrics = self._generic_epoch_end(training_step_outputs)
        self.log_metrics(metrics, split="train")

    def validation_step(self, batch, batch_idx):
        """Runs a prediction step for validation, logging the loss."""
        outputs = self._generic_step(batch, batch_idx)
        self.log(
            "loss/val",
            outputs["loss"],
            on_step=True,
            on_epoch=True,
            prog_bar=True,
            logger=True,
        )
        self.validation_step_outputs.append(outputs)
        return outputs

    def on_validation_epoch_end(self):
        if not self.trainer.sanity_checking:
            metrics = self._generic_epoch_end(self.validation_step_outputs)
            self.val_metrics = metrics  # cache for use in callback
            self.log_metrics(metrics, split="val")
            
    def test_step(self, batch, batch_idx):
        """Runs a predictionval_ step for testing, logging the loss."""
        outputs = self._generic_step(batch, batch_idx)
        self.test_step_outputs.append(outputs)
        return outputs

    def on_test_epoch_end(self):
        metrics = self._generic_epoch_end(self.test_step_outputs)
        self.test_metrics = metrics
        self.log_metrics(metrics, split="test")

    def log_metrics(self, metrics: typing.Dict, split: str):
        """Logs all metrics to logs and to tensorboard."""
        assert split in ["train", "val", "test"]

        # log our metrics to the logs directly
        log.info(f"{split.capitalize()} Accuracy: {metrics['accuracy']:.4f}")
        log.info(f"{split.capitalize()} Loss: {metrics['loss']:.4f}")

        # log metrics to tensorboard
        if split in ["train", "val"]:
            self.log(f"accuracy/{split}", metrics["accuracy"], on_epoch=True)
            self.log(f"loss/{split}", metrics["loss"], on_epoch=True)

        # log our metrics to the logs directly
        #metrics_summary = _summarize_metrics(
        #    metrics,
        #    self.class_names,
        #    split,
        #    self.current_epoch,
        #)
        #log.info(metrics_summary)

        # log metrics to tensorboard
        #if split in ["train", "val"]:
        #    self.log(f"precision/{split}", metrics["precision"], on_epoch=True)
        #    self.log(f"recall/{split}", metrics["recall"], on_epoch=True)
        #    self.log(f"f1_score/{split}", metrics["f1_score"], on_epoch=True)


    def predict_proba(self, features):
        """
        Returns the predicted probabilities for the given inputs.
        """
        logits = self.model(features)
        probabilities = torch.softmax(logits, dim=1)  # Apply softmax to get probabilities
        return probabilities.detach()
