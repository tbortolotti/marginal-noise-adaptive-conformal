import numpy as np
from sklearn.model_selection import train_test_split
import os
import hydra
import json
import torch
import pytorch_lightning as pl
from hydra.utils import instantiate
from hydra.core.hydra_config import HydraConfig
from hydra.core.global_hydra import GlobalHydra
from omegaconf import DictConfig, OmegaConf
from pathlib import Path

import torch.nn.functional as F
#from torchvision import transforms

from sklearn.metrics import accuracy_score, confusion_matrix

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm
import pdb

import sys, os
sys.path.append("..")
sys.path.append("../third_party")

from cln import data
from cln import contamination
from cln import estimation
from cln.utils import evaluate_predictions, estimate_rho
from cln.classification import MarginalLabelNoiseConformal

from third_party import arc
from third_party import bigearthnet

from third_party.bigearthnet.datamodules.bigearthnet_datamodule import BigEarthNetDataModule
from third_party.bigearthnet.models.bigearthnet_module import BigEarthNetModule

# Define default parameters
batch_size = 2000
epsilon_n_clean = 0.1
epsilon_n_corr = 0.1
estimate = "none"
seed = 1

# Parameters of the contamination process
#epsilon = 0.017
epsilon = 0.1
nu = 0.4
rho = np.array([0.15, 0.85])
contamination_model = "RRB"

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 9:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()

    batch_size = int(sys.argv[1])
    epsilon_n_clean = float(sys.argv[2])
    epsilon_n_corr = float(sys.argv[3])
    estimate = sys.argv[4]
    contamination_model = sys.argv[5]
    epsilon = float(sys.argv[6])
    nu = float(sys.argv[7])
    seed = int(sys.argv[8])


# Define other constant parameters
exp_num=201
data_name = "bigearthnet"
K = 5
epsilon_n = epsilon_n_clean + epsilon_n_corr
n_test = 500
num_exp = 5
allow_empty = True
epsilon_max = 0.1
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000

# Initialize noise contamination process
if contamination_model == "uniform":
    T = contamination.construct_T_matrix_simple(K, epsilon)
    M = contamination.convert_T_to_M(T,rho)
elif contamination_model == "block":
    T = contamination.construct_T_matrix_block(K, epsilon)
    M = contamination.convert_T_to_M(T,rho)
elif contamination_model == "RRB":
    T = contamination.construct_T_matrix_block_RR(K, epsilon, nu)
    M = contamination.convert_T_to_M(T,rho)
elif contamination_model == "random":
    T = contamination.construct_T_matrix_random(K, epsilon, random_state=seed)
    M = contamination.convert_T_to_M(T,rho)
else:
    print("Unknown contamination (M) model!")
    sys.stdout.flush()
    exit(-1)

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

# Load the pre-trained model
black_box = BigEarthNetModule(cfg)
mod_dir = cfg.out_directory.dir
mod_name = "trained_model.pth"
mod_path = os.path.join(mod_dir, mod_name)
black_box.load_state_dict(torch.load(mod_path))

# instantiate the datamodule
datamodule = BigEarthNetDataModule(
    cfg.datamodule.dataset_dir,
    cfg.datamodule.dataset_name,
    batch_size,
    cfg.datamodule.num_workers,
    transforms,
    label_mapping
)
datamodule.setup()

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'K':[K],
                       'n_cal':[n_cal], 'n_test':[n_test],
                       'epsilon_n_clean':[epsilon_n_clean], 'epsilon_n_corr':[epsilon_n_corr],
                       'estimate':[estimate], 'contamination':[contamination_model],
                       'epsilon':[epsilon], 'nu':[nu], 'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_n" + str(batch_size)
outfile_prefix += "_encl" + str(epsilon_n_clean) + "_enco" + str(epsilon_n_corr)
outfile_prefix += "_est" + estimate + "_" + "_cont" + contamination_model
outfile_prefix += "_eps" + str(epsilon) + "_nu" + str(nu) + "_" + str(seed)
print("Output file: {:s}.".format("results/"+outfile_prefix), end="\n")
sys.stdout.flush()

# Describe the experiment
def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    # Generate a large data set
    print("\nGenerating data...", end=' ')
    sys.stdout.flush()

    # Get reproducible random samples
    dataloader_train = datamodule.train_dataloader()
    dataloader_val = datamodule.val_dataloader()
    dataloader_test = datamodule.test_dataloader()

    batch = next(iter(dataloader_train))
    X_batch_train = batch['data']
    Y_batch_train = batch['labels']

    batch = next(iter(dataloader_val))
    X_batch_val = batch['data']
    Y_batch_val = batch['labels']

    batch = next(iter(dataloader_test))
    X_batch_test = batch['data']
    Y_batch_test = batch['labels']

    # Stack all features and labels together
    X_batch = torch.cat((X_batch_train, X_batch_val, X_batch_test), dim=0)
    Y_batch = torch.cat((Y_batch_train, Y_batch_val, Y_batch_test), dim=0)
    Y_batch = Y_batch.detach().numpy()
    print("Done.")
    sys.stdout.flush()

    ## TO DO ##
    #' Load Y_batch from file.... deve essere un numpy array
    # Y_batch = load from file i clean labels, oppure aggiungi una qualche corruption di qualche tipo
    #Y_batch = Yt_batch

    # Generate the contaminated labels
    print("Generating contaminated labels...", end=' ')
    sys.stdout.flush()
    contamination_process = contamination.LinearContaminationModel(T, random_state=random_state+3)
    Yt_batch = contamination_process.sample_labels(Y_batch)
    print("Done.")
    sys.stdout.flush()

    # Estimate the label proportions from the whole data set
    rho = estimate_rho(Y_batch, K)
    rho_tilde = estimate_rho(Yt_batch, K)

    # Separate the test set
    X, X_test, Y, Y_test, Yt, _ = train_test_split(X_batch, Y_batch, Yt_batch, test_size=n_test, random_state=random_state+2)

    # Estimate (if applicable) the label contamination model
    if estimate=="none":
        rho_hat = rho
        rho_tilde_hat = rho_tilde
        T_hat = contamination.construct_T_matrix_simple(K, epsilon)
        M_hat = contamination.convert_T_to_M(T_hat, rho_hat)
        epsilon_ci = None
        epsilon_hat = np.nan
    elif estimate=="rho":
        rho_tilde_hat = estimate_rho(Yt, K)
        T_hat = contamination.construct_T_matrix_simple(K, epsilon)  
        M_hat = contamination.convert_T_to_M(T_hat,rho)
        rho_hat = np.dot(M_hat.T, rho_tilde_hat)
        epsilon_ci = None
        epsilon_hat = np.nan
    elif estimate=="rho-epsilon-point":
        # Hold-out some data to estimate the contamination model
        X, X_estim, Y, Y_estim, Yt, Yt_estim = train_test_split(X, Y, Yt, test_size=epsilon_n, random_state=random_state+3)

        # Keep some hold-out data clean
        X_estim_clean, X_estim_corr, Y_estim_clean, _, _, Yt_estim_corr = train_test_split(X_estim, Y_estim, Yt_estim,
                                                                                           test_size=epsilon_n_corr/epsilon_n, random_state=random_state+4)

        rho_tilde_hat = estimate_rho(Yt, K)
        epsilon_hat, _, _, _, _ = estimation.fit_contamination_model_RR(X_estim_clean, X_estim_corr,
                                                                        Y_estim_clean, Yt_estim_corr, black_box,
                                                                        K, 0.01, pre_trained=True,
                                                                        random_state=random_state+6)
        T_hat = contamination.construct_T_matrix_simple(K, epsilon_hat)
        rho_hat = np.dot(np.linalg.inv(T_hat), rho_tilde_hat)
        M_hat = contamination.convert_T_to_M(T_hat,rho_hat)
        epsilon_ci = None

    else:
        print("Unknown estimation option!")
        sys.stdout.flush()
        exit(-1)


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

                "Adaptive": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                epsilon=epsilon, T=T_hat, M=M_hat, rho_tilde=rho_tilde_hat,
                                                                allow_empty=allow_empty, method="old", optimistic=False,
                                                                verbose=False, pre_trained=True, random_state=random_state),

                "Adaptive+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box, K, alpha, n_cal=-1,
                                                                epsilon=epsilon, T=T_hat, M=M_hat, rho_tilde=rho_tilde_hat,
                                                                allow_empty=allow_empty, method="old", optimistic=True,
                                                                verbose=False, pre_trained=True, random_state=random_state),

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
                                                                   pre_trained=True, random_state=random_state)

            }

            # Initialize an empty list to store the evaluation results
            res_list = []

            # Loop through the methods, apply them, and evaluate the results
            for method_name, method_func in methods.items():
                print(f"Applying {method_name} method...", end=' ')
                sys.stdout.flush()

                # Initialize and apply the method
                method = method_func()
                predictions = method.predict(X_test)

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
