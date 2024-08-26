import os
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from torchvision.models import resnet50
from BHN_dataset import BigEarthNet

import pdb
import pickle

tar_path = os.path.expanduser('~/data/BigEarthNet/BigEarthNet-S2-v1.0.tar.gz')
df = pd.read_csv('BHN_labels.csv')

df_train = df[df['split']=='train'].copy()
train_dataset = BigEarthNet(df_train, tar_path, transform=True)
train_size = len(train_dataset)

df_test = df[df['split']=='test'].copy()
test_dataset = BigEarthNet(df_test, tar_path, transform=True)
test_size = len(test_dataset)

K = df['v1-labels'].nunique()
epochs = 20
seed = 1
batch_size = 64

torch.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
np.random.seed(seed)

#train_size = 10000
#test_size = 20000
#train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

#pdb.set_trace()

# 3. Training the net (ResNet-50 modified for 12 channels)
model = resnet50(weights=None)
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

# Save the model's state dictionary
model_save_path = 'resnet50_trained.pth'
torch.save(model.state_dict(), model_save_path)

# Store relevant parameters
info_save = {
    'num_epochs': epochs,
    'batch_size': batch_size,
    'K': K,
    'seed': seed,
    'accuracy': accuracy
}

filename = f'training_parameters_{epochs}_epochs.pkl'
# Save the dictionary to a file
with open('filename', 'wb') as f:
    pickle.dump(info_save, f)



