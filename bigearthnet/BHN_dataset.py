import tarfile
import random

import os
import random
from PIL import Image
import numpy as np
import pandas as pd
import cv2
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import torchvision.transforms as transforms
from torchvision.models import resnet50

bands_stats = {
    'mean': {
        'B01': 340.76769064,
        'B02': 429.9430203,
        'B03': 614.21682446,
        'B04': 590.23569706,
        'B05': 950.68368468,
        'B06': 1792.46290469,
        'B07': 2075.46795189,
        'B08': 2218.94553375,
        'B8A': 2266.46036911,
        'B09': 2246.0605464,
        'B11': 1594.42694882,
        'B12': 1009.32729131
    },
    'std': {
        'B01': 554.81258967,
        'B02': 572.41639287,
        'B03': 582.87945694,
        'B04': 675.88746967,
        'B05': 729.89827633,
        'B06': 1096.01480586,
        'B07': 1273.45393088,
        'B08': 1365.45589904,
        'B8A': 1356.13789355,
        'B09': 1302.3292881,
        'B11': 1079.19066363,
        'B12': 818.86747235
    }
}

labels = {
    'Urban fabric',
    'Mixed',
    'Non-Urban'
}

def normalize(img, mean, std):
    min_value = mean - 2 * std
    max_value = mean + 2 * std
    img = (img - min_value) / (max_value - min_value) * 255.0
    img = np.clip(img, 0, 255).astype(np.uint8)
    return img

class BigEarthNet(Dataset):
    def __init__(self, df, tar_path, transform_flag=False):
        self.df = df
        self.tar_path = tar_path
        self.transform_flag = transform_flag
        self.bands = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B11', 'B12']
    
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        patch_id = self.df.iloc[idx]['patch_id']
        with tarfile.open(self.tar_path, 'r:gz') as tar:
            # List of files associated to the patch_id
            band_files = [f for f in tar.getnames() if patch_id in f and f.endswith('.tif')]
            # Sort files according to band
            band_files.sort()

            # Read all channels of the image
            bands = []
            for b, band_file in enumerate(band_files):
                band = tar.extractfile(band_file)
                image = Image.open(band)
                image_arr = np.array(image)
                image_arr = normalize(image_arr, mean=bands_stats['mean'][self.bands[b]], std=bands_stats['std'][self.bands[b]])
                image_arr = cv2.resize(image_arr, dsize=(128, 128), interpolation=cv2.INTER_CUBIC)
                bands.append(np.array(image_arr))

            # Combine the channels in a single image
            image = np.stack(bands, axis=-1)

        if self.transform_flag:
            transform = transforms.Compose([
                transforms.ToPILImage(),
                #transforms.Resize((224, 224)),
                transforms.Resize((128, 128)),
                transforms.CenterCrop(112),
                transforms.ToTensor()])
            image = transform(image)
        
        label = self.df.loc[self.df['patch_id']==patch_id, 'v1-label'].values[0]

        return image, label