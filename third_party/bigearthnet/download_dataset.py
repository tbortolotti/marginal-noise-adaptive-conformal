import pathlib
from datamodules.bigearthnet_datamodule import download_data

# Define dataset directory and name
dataset_dir = "../datasets"  # Adjust this to your desired path
dataset_name = "bigearthnet-full"

# Download the dataset
download_data(dataset_dir, dataset_name)
