import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt

import torch
import torch.nn.functional as F
from torchvision import transforms as T
from torchvision.datasets import CIFAR10
import torch.nn as nn
from torchvision.models import resnet18 as imagenet_resnet18, ResNet18_Weights

import sys, os

import_path = "/home1/tb_214/code/PyTorch_CIFAR10"
sys.path.append(import_path)

from cifar10_models.resnet import resnet18 as cifar_resnet18

import pickle
def unpickle(file):
    import pickle
    with open(file, 'rb') as fo:
        dict = pickle.load(fo, encoding='latin1')
    return dict

import pdb


def draw_images(images, labels, rows=5, columns = 5):
    images = images.detach().numpy()
    images = images.reshape(len(images),3,32,32).transpose(0,2,3,1)
    assert len(images) >= rows*columns
    fig=plt.figure(figsize=(10, 10))
    # visualize these random images
    for i in range(1, columns*rows +1):
        fig.add_subplot(rows, columns, i)
        plt.imshow(images[i-1])
        plt.xticks([])
        plt.yticks([])
        plt.title("{}".format(labels[i-1]))
    plt.show()


class Cifar10DataSet:

    def __init__(self, data_dir, noisy_data_dir, normalize=True, random_state=2023):

        self.rng = np.random.default_rng(seed=random_state)

        MEAN = [0.4914, 0.4822, 0.4465]
        STD = [0.2471, 0.2435, 0.2616]

        if normalize:
            transform = T.Compose(
                [
                    T.ToTensor(),
                    T.Normalize(MEAN, STD),
                ]
            )
        else:
            transform = T.Compose(
                [
                    T.ToTensor(),
                ]
            )

        self.cifar10 = CIFAR10(root=data_dir, train=False, transform=transform)

        meta_file = data_dir + "/cifar-10-batches-py/batches.meta"
        meta_data = unpickle(meta_file)
        self.label_names = np.array(meta_data['label_names'])

        # Load noisy labels
        noisy_label_file = noisy_data_dir + "/cifar10h-counts.npy"
        self.noisy_label_data = np.load(noisy_label_file)

    def generate_noisy_labels(self, index):
        weights = self.noisy_label_data[index]
        n,K = weights.shape
        weights = weights / np.sum(weights,1).reshape(n,1)
        noisy_labels = [self.rng.choice(K, replace=True, p=weights[i]) for i in range(n)]
        noisy_labels_name = self.label_names[noisy_labels]
        return noisy_labels[0], list(noisy_labels_name)[0]

    def __getitem__(self, index):
        data, label = self.cifar10[index]
        label_name = self.label_names[label]
        noisy_label, noisy_label_name = self.generate_noisy_labels([index])
        return data, label, label_name, noisy_label, noisy_label_name, index

    def __len__(self):
        return len(self.cifar10)


class ResNet18:
    def __init__(self):
        self.black_box = cifar_resnet18(pretrained=True)
        self.black_box.eval()

    def predict(self, X):
        pi_hat = self.predict_proba(X)
        Y_hat = np.argmax(pi_hat, 1)
        return Y_hat

    def predict_proba(self, X):
        Z = self.black_box(X)
        pi_hat = torch.softmax(Z, dim=1).detach().numpy()
        return pi_hat
    
class ImageNetResNet18Features:
    def __init__(self, device=None):
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")

        weights = ResNet18_Weights.DEFAULT
        model = imagenet_resnet18(weights=weights)

        # Remove final classification layer
        self.backbone = nn.Sequential(*list(model.children())[:-1])
        self.backbone.to(self.device)
        self.backbone.eval()

        # Freeze parameters
        for p in self.backbone.parameters():
            p.requires_grad = False

    def transform(self, X):
        X = X.to(self.device)
        with torch.no_grad():
            feats = self.backbone(X)              # (N, 512, 1, 1)
            feats = feats.view(feats.size(0), -1) # (N, 512)
        return feats.cpu()

class Cifar10ImageNetDataSet:
    def __init__(self, data_dir, noisy_data_dir, random_state=2026):
        self.rng = np.random.default_rng(seed=random_state)

        # ImageNet preprocessing
        self.transform = T.Compose([
            T.Resize((224, 224)),
            T.ToTensor(),
            T.Normalize(
                mean=[0.485, 0.456, 0.406],
                std=[0.229, 0.224, 0.225]
            ),
        ])

        self.cifar10 = CIFAR10(root=data_dir, train=False, transform=self.transform)

        meta_file = data_dir + "/cifar-10-batches-py/batches.meta"
        meta_data = unpickle(meta_file)
        self.label_names = np.array(meta_data["label_names"])

        noisy_label_file = noisy_data_dir + "/cifar10h-counts.npy"
        self.noisy_label_data = np.load(noisy_label_file)

    def generate_noisy_labels(self, index):
        weights = self.noisy_label_data[index]
        n, K = weights.shape
        weights = weights / np.sum(weights, 1).reshape(n, 1)
        noisy_labels = [self.rng.choice(K, replace=True, p=weights[i]) for i in range(n)]
        noisy_labels_name = self.label_names[noisy_labels]
        return noisy_labels[0], list(noisy_labels_name)[0]

    def __getitem__(self, index):
        data, label = self.cifar10[index]
        label_name = self.label_names[label]
        noisy_label, noisy_label_name = self.generate_noisy_labels([index])
        return data, label, label_name, noisy_label, noisy_label_name, index

    def __len__(self):
        return len(self.cifar10)
