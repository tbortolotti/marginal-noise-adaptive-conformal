import numpy as np
from sklearn.model_selection import train_test_split
import torch
import torch.nn.functional as F
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
from cln.T_estimation_NN import NoisyLabelNet, train_alternate
from cln.T_estimation import TMatrixEstimation
from cln.utils import evaluate_predictions, estimate_rho
from cln.classification import MarginalLabelNoiseConformal
from third_party import arc

from data_torch import Cifar10DataSet, CifarResNet18Features, ResNet18
from torch.utils.data import DataLoader
import gc
import copy


# Define default parameters
exp_num = 911
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n_train = 2000
n_clean = 500
n_cal = 1000
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 10:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()
    exp_num = int(sys.argv[1])
    epsilon = float(sys.argv[2])
    contamination_model = sys.argv[4]
    n_train = int(sys.argv[5])
    n_clean = int(sys.argv[6])
    n_cal = int(sys.argv[7])
    contamination_exp_flag = sys.argv[8].lower() == "true"
    seed = int(sys.argv[9])

# Define other constant parameters
data_name = "cifar10"
K = 10
num_exp = 5
allow_empty = True
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000
nu = 0.2
n_test = 500
batch_size = n_train + n_clean + n_cal + n_test
epsilon_init = 0

# Set default directories
data_dir = "/home1/tb_214/data/cifar10"
noisy_data_dir = "/home1/tb_214/data/cifar10h"

print(f"Data Directory: {data_dir}")
print(f"Noisy Data Directory: {noisy_data_dir}")

device = "cuda" if torch.cuda.is_available() else "cpu"

dataset = Cifar10DataSet(data_dir=data_dir, noisy_data_dir=noisy_data_dir, random_state=2026)

#print(f"\nOverall numerosity of dataset {len(dataset)}")
#sys.stdout.flush()

loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=1)
feature_extractor = CifarResNet18Features()

# Initialize noise contamination process
if contamination_model == "uniform":
    T = contamination.construct_T_matrix_simple(K, epsilon)
elif contamination_model == "block":
    T = contamination.construct_T_matrix_block(K, epsilon)
elif contamination_model == "RRB":
    T = contamination.construct_T_matrix_block_RR(K, epsilon, nu)

# Initialize black-box model
black_box = ResNet18()

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'K':[K],
                       'n_train':[n_train], 'n_clean':[n_clean], 'n_cal':[n_cal],
                       'epsilon':[epsilon], 'contamination':[contamination_model],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_n" + str(batch_size)
outfile_prefix += "_eps" + str(epsilon) + "_" + contamination_model
outfile_prefix += "_nt" + str(n_train) + "_ncl" + str(n_clean) +"_nc" + str(n_cal) + "_seed" + str(seed)
print("Output file: {:s}.".format("results/"+outfile_prefix), end="\n")
sys.stdout.flush()

# Describe the experiment
def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    # Generate dataset for calibration and test
    # Generate a large data set
    print("\nGenerating data...", end=' ')
    sys.stdout.flush()
    if contamination_model == "true":
        X_all, _, Y_all, _, Yt_all, _, _ = next(iter(loader))
        Y_all = Y_all.detach().numpy()
        Yt_all = Yt_all.detach().numpy()
        print("Done.")
        sys.stdout.flush()

        # Separate the test set
        X, X_test, Y, Y_test, Yt, Yt_test = train_test_split(X_all, Y_all, Yt_all, test_size=n_test, random_state=random_state+1)
        del X_all, Y_all, Yt_all, Yt_test
    else:
        X_all, _, Y_all, _, _, _, _ = next(iter(loader))
        Y_all = Y_all.detach().numpy()
        print("Done.")
        sys.stdout.flush()

        # Separate the test set
        X, X_test, Y, Y_test = train_test_split(X_all, Y_all, test_size=n_test, random_state=random_state+1)
        del X_all, Y_all

        # Generate the contaminated labels
        print("Generating contaminated labels...", end=' ')
        sys.stdout.flush()
        contamination_process = contamination.LinearContaminationModel(T, random_state=random_state+2)
        Yt = contamination_process.sample_labels(Y)
        print("Done.")
        sys.stdout.flush()


    # Estimate the label proportions from the whole data set
    print("Estimating label proportions...", end=' ')
    sys.stdout.flush()
    rho_tilde_hat = estimate_rho(Yt, K)
    #print(rho_tilde_hat)
    print("Done.")
    sys.stdout.flush()

    # Separate data into training and calibration
    X_train, X_cal, Y_train, Y_cal, Yt_train, Yt_cal = train_test_split(X, Y, Yt, test_size=n_cal, random_state=random_state+3)
    del X, Y, Yt

    print("Generating clean dataset...", end=' ')
    sys.stdout.flush()
    # Identify the central observations in the training set to build the clean dataset
    conf_scores = black_box.predict_proba(X_train).max(axis=1)
    top_indices = np.argsort(conf_scores)[-n_clean:]

    X_clean = X_train[top_indices]
    Y_clean = Y_train[top_indices]
    Yt_clean = Yt_train[top_indices]
    I = np.zeros(len(Y_train))
    I[top_indices] = 1
    Y_obs = np.where(I == 1, Y_train, Yt_train)
    print("Done.")
    sys.stdout.flush()

    # Extract features
    print("Extract features for T estimation...", end=' ')
    sys.stdout.flush()

    X_cifar_feat = feature_extractor.transform(X_train)
    num_var_cifar = X_cifar_feat.shape[1]
    X_cifar_feat_torch = X_cifar_feat.to(device)
    Y_obs_torch = torch.tensor(Y_obs, dtype=torch.long)
    I_torch = torch.tensor(I, dtype=torch.long)

    del X_train, X_cifar_feat

    # Use extracted features to estimate the contamination
    if not contamination_exp_flag:
        #____________________________________________________________________
        ## Estimate T using the clean/noisy correspondence
        T_method = TMatrixEstimation(Y_clean, Yt_clean, K, estimation_method="empirical_parametricRR")
        T_hat_clean = T_method.get_estimate()

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[64], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features and SLL...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN_cifar_sll = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar_sll = T_hat_NN_cifar_sll.detach().numpy()
        print("Done.")
        sys.stdout.flush()

    else:
        #____________________________________________________________________
        ## Estimate T using the clean/noisy correspondence
        T_method = TMatrixEstimation(Y_clean, Yt_clean, K, estimation_method="empirical")
        T_hat_clean = T_method.get_estimate()

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[64], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features and SLL...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN_cifar_sll = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar_sll = T_hat_NN_cifar_sll.detach().numpy()
        print("Done.")
        sys.stdout.flush()

    # Force garbage collection
    gc.collect()
    torch.cuda.empty_cache()

    """
    # Print memory state
    print(f"RAM used: {torch.cuda.memory_allocated()/1e9:.2f} GB")
    print(f"RAM reserved: {torch.cuda.memory_reserved()/1e9:.2f} GB")

    # Print all live tensors
    for obj in gc.get_objects():
        try:
            if torch.is_tensor(obj) or (hasattr(obj, 'data') and torch.is_tensor(obj.data)):
                print(f"Tensor: {type(obj).__name__}, size={obj.size()}, device={obj.device}")
        except:
            pass
    """

    alpha = 0.1
    guarantee = 'marginal'

    res = pd.DataFrame({})

    print("\nSeeking {:s} coverage at level {:.2f}.".format(guarantee, 1-alpha))

    # Define a dictionary of methods with their names and corresponding initialization parameters
    methods = {
        "Standard": lambda: arc.methods.SplitConformal(X_cal, Yt_cal, black_box, K, alpha, n_cal=-1,
                                                       pre_trained=True, random_state=random_state),

        "Standard using clean": lambda: arc.methods.SplitConformal(X_clean, Y_clean, black_box, K, alpha, n_cal=-1,
                                                                pre_trained=True, random_state=random_state),

        "Adaptive optimized+": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ clean": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_clean, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ NN": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_NN_cifar, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ NN SLL": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_NN_cifar_sll, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
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
