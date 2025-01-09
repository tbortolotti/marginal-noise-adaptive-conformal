import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import os
import hydra
import json
import torch
import pytorch_lightning as pl
from hydra.utils import instantiate
from hydra.core.global_hydra import GlobalHydra
from pathlib import Path

import torch.nn.functional as F

from sklearn.metrics import accuracy_score, confusion_matrix

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
#from matplotlib import pyplot as plt
#from tqdm import tqdm
#import pdb

import sys, os
sys.path.append("..")
sys.path.append("../third_party")

from cln import contamination
from cln.utils import evaluate_predictions
from cln.classification import MarginalLabelNoiseConformal
from cln.classification_label_conditional import LabelNoiseConformal

from third_party import arc
from third_party import bigearthnet
from third_party.bigearthnet.datamodules.bigearthnet_datamodule import BigEarthNetDataModule
from third_party.bigearthnet.models.bigearthnet_module import BigEarthNetModule

# Define default parameters
batch_size = 2000
estimate = "none"
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 5:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()
    exp_num=int(sys.argv[1])
    batch_size = int(sys.argv[2])
    estimate = sys.argv[3]
    seed = int(sys.argv[4])

# Define other constant parameters
data_name = "bigearthnet"
n_test = 500
num_exp = 5
allow_empty = True
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000

# Oracle parameters
rho = [0.114, 0.032, 0.026, 0.138, 0.001, 0.689]
rho_tilde = [0.113, 0.031, 0.025, 0.137, 0.016, 0.678]
epsilon = 0.016
K = 6
T = np.array([
    [0.986, 0, 0, 0, 0, 0],
    [0, 0.988, 0, 0, 0, 0],
    [0, 0, 0.956, 0, 0, 0],
    [0, 0, 0, 0.994, 0, 0.001],
    [0.013, 0.008, 0.042, 0.006, 1, 0.017],
    [0.001, 0.004, 0.002, 0, 0, 0.982]
])

# Pre-process parameters
n_cal = batch_size - n_test

with open('../third_party/bigearthnet/data/label_mapping.json', 'r') as f:
    label_mapping = json.load(f)

# Clear any previous Hydra instances if they exist
if GlobalHydra.instance().is_initialized():
    GlobalHydra.instance().clear()
hydra.initialize(config_path="../third_party/bigearthnet/configs", version_base="1.2")

# fetch the transforms used in the model
cfg = hydra.compose(config_name="config")
transforms = instantiate(cfg.transforms.obj)

# Load clean labels
dataset_name = cfg.datamodule.dataset_name
file_path_train = os.path.join(
    '../third_party/bigearthnet/data',
    f'train_{dataset_name}.csv'
)
v1v2_corresp_train = pd.read_csv(file_path_train, header=0)

"""
file_path_val = os.path.join(
    '../third_party/bigearthnet/data',
    f'val_{dataset_name}.csv'
)
v1v2_corresp_val = pd.read_csv(file_path_val, header=0)

file_path_test = os.path.join(
    '../third_party/bigearthnet/data',
    f'test_{dataset_name}.csv'
)
v1v2_corresp_test = pd.read_csv(file_path_test, header=0)
"""

# Load the pre-trained model
black_box = BigEarthNetModule(cfg)
mod_dir = cfg.out_directory.dir
mod_name = "trained_model.pth"
mod_path = os.path.join(mod_dir, mod_name)
black_box.load_state_dict(torch.load(mod_path))

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'K':[K],
                       'n_cal':[n_cal], 'n_test':[n_test],
                       'estimate':[estimate],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_n" + str(batch_size)
outfile_prefix += "_est" + estimate + "_" + str(seed)
print("Output file: {:s}.".format("results/"+outfile_prefix), end="\n")
sys.stdout.flush()

# Describe the experiment
def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    # Generate a large data set
    print("\nGenerating data...", end=' ')
    sys.stdout.flush()

    # instantiate the datamodule
    datamodule = BigEarthNetDataModule(
    cfg.datamodule.dataset_dir,
    cfg.datamodule.dataset_name,
    batch_size,
    cfg.datamodule.num_workers,
    transforms,
    label_mapping,
    int(seed*random_state),
    )
    datamodule.setup()

    # Get reproducible random samples
    dataloader_train = datamodule.train_dataloader()

    batch = next(iter(dataloader_train))
    X_batch_train = batch['data']
    Yt_batch_train = batch['labels']
    generator = torch.Generator().manual_seed(datamodule.random_seed)
    indices_df = torch.randperm(len(datamodule.train_dataset), generator=generator).tolist()
    shuffled_csv_df = v1v2_corresp_train.iloc[indices_df].reset_index(drop=True)
    batch_df = shuffled_csv_df.iloc[0 : int(batch_size)]
    Y_batch_train = batch_df['v2-labels-grouped'].to_numpy()
    valid_indices = torch.tensor(~np.isnan(Y_batch_train), dtype=torch.bool)
    X_batch_train = X_batch_train[valid_indices,:,:,:]
    Yt_batch_train = Yt_batch_train[valid_indices]
    Y_batch_train = Y_batch_train[valid_indices].astype(int)

    # Stack all features and labels together
    X_batch = X_batch_train
    Yt_batch = Yt_batch_train
    Yt_batch = Yt_batch.detach().numpy()
    Y_batch = Y_batch_train
    print(f"Done. The dimension of the current batch is: {len(Yt_batch)}")
    sys.stdout.flush()

    # Separate the test set
    X, X_test, Y, Y_test, Yt, _ = train_test_split(X_batch, Y_batch, Yt_batch, test_size=n_test, random_state=random_state+2)

    # Plug-in estimate of the label contamination model
    rho_tilde_hat = rho_tilde
    rho_hat = rho
    T_hat = T
    M_hat = contamination.convert_T_to_M(T_hat, rho_hat)

    res = pd.DataFrame({})
    for alpha in [0.1]:
        for guarantee in ['marginal']:

            print("\nSeeking {:s} coverage at level {:.2f}.".format(guarantee, 1-alpha))

            label_conditional = False
            alpha_theory = alpha * (1 - epsilon * (1-1/K))

            # Define a dictionary of methods with their names and corresponding initialization parameters
            methods = {
                "Standard": lambda: arc.methods.SplitConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                               label_conditional=label_conditional, allow_empty=allow_empty,
                                                               pre_trained=True, random_state=random_state),
                
                "Standard (theory)": lambda: arc.methods.SplitConformal(X, Yt, black_box, K, alpha_theory, n_cal=-1,
                                                               label_conditional=label_conditional, allow_empty=allow_empty,
                                                               pre_trained=True, random_state=random_state),

                "Adaptive optimized": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                          epsilon=epsilon, T=T_hat, rho_tilde=rho_tilde_hat,
                                                                          allow_empty=allow_empty, method="improved",
                                                                          optimized=True, optimistic=False, verbose=False,
                                                                          pre_trained=True, random_state=random_state),

                "Adaptive optimized+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                          epsilon=epsilon, T=T_hat, rho_tilde=rho_tilde_hat,
                                                                          allow_empty=allow_empty, method="improved",
                                                                          optimized=True, optimistic=True, verbose=False,
                                                                          pre_trained=True, random_state=random_state),

                "Adaptive simplified": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                           epsilon=epsilon, T=T_hat, rho_tilde=rho_tilde_hat,
                                                                           allow_empty=allow_empty, method="improved",
                                                                           optimized=False, optimistic=False, verbose=False,
                                                                           pre_trained=True, random_state=random_state),

                "Adaptive simplified+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                           epsilon=epsilon, T=T_hat, rho_tilde=rho_tilde_hat,
                                                                           allow_empty=allow_empty, method="improved",
                                                                           optimized=False, optimistic=True, verbose=False,
                                                                           pre_trained=True, random_state=random_state),

                "Asymptotic": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                  epsilon=epsilon, asymptotic_h_start=asymptotic_h_start,
                                                                  asymptotic_MC_samples=asymptotic_MC_samples, T=T_hat,
                                                                  rho_tilde=rho_tilde_hat, allow_empty=allow_empty,
                                                                  method="asymptotic", optimistic=False, verbose=False,
                                                                  pre_trained=True, random_state=random_state),

                "Asymptotic+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                   epsilon=epsilon, asymptotic_h_start=asymptotic_h_start,
                                                                   asymptotic_MC_samples=asymptotic_MC_samples, T=T_hat,
                                                                   rho_tilde=rho_tilde_hat, allow_empty=allow_empty,
                                                                   method="asymptotic", optimistic=True, verbose=False,
                                                                   pre_trained=True, random_state=random_state),

                "Label conditional": lambda: LabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                 rho_tilde=rho_tilde_hat, M=M_hat,
                                                                 calibration_conditional=False, gamma=None,
                                                                 optimistic=False, allow_empty=allow_empty, verbose=False, pre_trained=True, random_state=random_state),
                
                "Label conditional+": lambda: LabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                  rho_tilde=rho_tilde_hat, M=M_hat,
                                                                  calibration_conditional=False, gamma=None,
                                                                  optimistic=True, allow_empty=allow_empty, verbose=False, pre_trained=True, random_state=random_state)

            }

            # Initialize an empty list to store the evaluation results
            res_list = []

            # Loop through the methods, apply them, and evaluate the results
            for method_name, method_func in methods.items():
                print(f"Applying {method_name} method...", end=' ')
                sys.stdout.flush()

                # Initialize and apply the method
                method = method_func()
                predictions = method.predict(X_test, random_state=2023)

                print("Done.")
                sys.stdout.flush()

                # Evaluate the method
                res_new = evaluate_predictions(predictions, X_test, Y_test, K, verbose=False)
                res_new['Method'] = method_name
                res_new['Guarantee'] = guarantee
                res_new['Alpha'] = alpha
                res_new['random_state'] = random_state

                # Append the result to the results list
                res_list.append(res_new)


            # Combine all results into a single DataFrame
            res = pd.concat(res_list)
                
    print(res)

    return res

# Run all experiments
results = pd.DataFrame({})
for batch in np.arange(1,num_exp+1):
    res = run_experiment(1000*seed+batch-1000)
    results = pd.concat([results, res])

    # Save results
    outfile = "results/" + outfile_prefix + ".txt"
    results_out = pd.concat([header,results], axis=1)
    results_out.to_csv(outfile, index=False, float_format="%.5f")

print("\nPreview of results:")
print(results)
sys.stdout.flush()

print("\nSummary of results:")
summary = results.groupby(['Alpha', 'Guarantee', 'Method', 'Label']).agg(['mean','std']).reset_index()
print(summary)
sys.stdout.flush()


print("\nFinished.\nResults written to {:s}\n".format(outfile))
sys.stdout.flush()
