import tarfile
import random

import os
import random
from PIL import Image
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import torchvision.transforms as transforms
from torchvision.models import resnet50
from BHN_dataset import BigEarthNet


tar_path = "/data/BigEarthNet/BigEarthNet-S2-v1.0.tar.gz"
df = pd.read_csv('BHV_labels.csv')

dataset = BigEarthNet(df, tar_path, transform=True)

"""
K = 3
epochs = 10

train_size = 10000
test_size = 10000
train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size], random_state=42)
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

# 3. Training the net (ResNet-50 modified for 12 channels)
model = resnet50(pretrained=False)
model.conv1 = nn.Conv2d(
    in_channels=12,
    out_channels=64,
    kernel_size=(7, 7),
    stride=(2, 2),
    padding=(3, 3),
    bias=False
)

model.fc = nn.Linear(model.fc.in_features, K)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)

model.fc.weight.data.normal_(mean=0.0,std=0.01)
model.fc.bias.data.zero_()

criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

for epoch in range(epochs):
    model.train()
    for images, labels in train_loader:
        images = images.to(device)
        labels = labels.to(device)
        
        optimizer.zero_grad()
        outputs = model(images)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

# Prediction error
model.eval()
correct = 0
total = 0
with torch.no_grad():
    for images, labels in test_loader:
        images = images.to(device)
        labels = labels.to(device)
        outputs = model(images)
        _, predicted = torch.max(outputs.data, 1)
        total += labels.size(0)
        correct += (predicted == labels).sum().item()

accuracy = correct / total

"""
