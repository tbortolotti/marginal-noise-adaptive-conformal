import numpy as np
import torch
#from torchvision import transforms

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm
import pdb

import sys, os
sys.path.append("..")
sys.path.append("../third_party")

from cln import contamination
from third_party import arc

#from data_torch import Cifar10DataSet, ResNet18
from data_torch import Cifar10ImageNetDataSet, ImageNetResNet18Features
from torch.utils.data import DataLoader

# Define default parameters
exp_num = 1000
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n_train1 = 1000
n_train2 = 1000
seed = 1

# Define other constant parameters
batch_size = n_train1 + n_train2
data_name = "cifar10"
K = 10
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000

# Set default directories
data_dir = "/home1/tb_214/data/cifar10"
noisy_data_dir = "/home1/tb_214/data/cifar10h"

print(f"Data Directory: {data_dir}")
print(f"Noisy Data Directory: {noisy_data_dir}")

dataset = Cifar10ImageNetDataSet(data_dir=data_dir, noisy_data_dir=noisy_data_dir, random_state=2026)
loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=1)

# Initialize the black-box model
#black_box_RN18 = ResNet18()
black_box_SVC = arc.black_boxes.SVC(clip_proba_factor = 1e-5)

# Initialize noise contamination process
T = contamination.construct_T_matrix_simple(K, epsilon)

print("\nRunning experiment in batch 1...")
sys.stdout.flush()

# Generate a large data set
print("\nGenerating data...", end=' ')
sys.stdout.flush()

feature_extractor = ImageNetResNet18Features()

X, Y, Y_lab, Yt, Yt_lab, idx_batch = next(iter(loader))
X_features = feature_extractor.transform(X)
X_features = X_features.numpy()

print(X_features.shape, len(Y), len(Yt))
sys.stdout.flush()

print("\nTraining the black box...", end=' ')
sys.stdout.flush()
black_box_SVC.fit(X_features, Yt)

p_hat = black_box_SVC.predict_proba(X_features)
print(p_hat.shape)
sys.stdout.flush()

print("\nNo error.... nice", end=' ')
sys.stdout.flush()

