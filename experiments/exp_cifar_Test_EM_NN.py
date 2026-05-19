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
from cln.T_estimation import evaluate_estimate
from cln.T_estimation_NN import NoisyLabelNet, train_alternate

#from data_torch import ResNet18
#from data_torch import Cifar10DataSet, ImageNetResNet18Features, CifarResNet18Features, ResNet18
from data_torch import Cifar10DataSet, CifarResNet18Features, ResNet18
from torch.utils.data import DataLoader
import gc

# Define default parameters
exp_num = 821
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n = 1000
n_clean = 500
pi_clean = 0
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
    nu = float(sys.argv[3])
    contamination_model = sys.argv[4]
    n = int(sys.argv[5])
    n_clean = int(sys.argv[6])
    pi_clean = float(sys.argv[7])
    contamination_exp_flag = sys.argv[8].lower() == "true"
    seed = int(sys.argv[9])


# Define other constant parameters
data_name = "cifar10"
K = 10
num_exp = 5
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000
epsilon_init = 0

if n_clean == 0:
    n_clean = int(np.round(pi_clean * n))
    n_noisy = n - n_clean
else:
    n_noisy = n
    n = n_noisy + n_clean

batch_size = n

# Set default directories
data_dir = "/home1/tb_214/data/cifar10"
noisy_data_dir = "/home1/tb_214/data/cifar10h"

print(f"Data Directory: {data_dir}")
print(f"Noisy Data Directory: {noisy_data_dir}")

device = "cuda" if torch.cuda.is_available() else "cpu"

dataset = Cifar10DataSet(data_dir=data_dir, noisy_data_dir=noisy_data_dir, random_state=2026)
loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=1)
#feature_extractor = ImageNetResNet18Features()
cifar_feature_extractor = CifarResNet18Features()

# Initialize the black-box model
black_box = ResNet18()

# Initialize noise contamination process
if contamination_model == "uniform":
    T = contamination.construct_T_matrix_simple(K, epsilon)
elif contamination_model == "block":
    T = contamination.construct_T_matrix_block(K, epsilon)
elif contamination_model == "RRB":
    T = contamination.construct_T_matrix_block_RR(K, epsilon, nu)
elif contamination_model == "random":
    T = contamination.construct_T_matrix_random(K, epsilon, random_state=seed)
else:
    print("Unknown contamination model!")
    sys.stdout.flush()
    exit(-1)

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'K':[K],
                       'n':[n], 'n_noisy':[n_noisy],
                       'n_clean':[n_clean], 'pi_clean':[pi_clean],
                       'epsilon':[epsilon], 'nu':[nu], 'contamination':[contamination_model],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_n" + str(batch_size)
outfile_prefix += "_eps" + str(epsilon) + "_nu" + str(nu) + "_" + contamination_model
outfile_prefix += "_n" + str(n) + "_ncl" + str(n_clean) + "_picl" + str(pi_clean) + "_seed" + str(seed)
print("Output file: {:s}.".format("results/"+outfile_prefix), end="\n")
sys.stdout.flush()

# Describe the experiment
def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    # Generate a large data set
    print("\nGenerating data...", end=' ')
    sys.stdout.flush()

    X, _, Y, _, _, _, _ = next(iter(loader))
    #X_features = feature_extractor.transform(X)
    #X_features = X_features.numpy()
    #num_var = X_features.shape[1]

    X_cifar_feat = cifar_feature_extractor.transform(X)
    num_var_cifar = X_cifar_feat.shape[1]
    
    Y = Y.detach().numpy()

    # Generate the contaminated labels
    print("Generating contaminated labels...", end=' ')
    sys.stdout.flush()
    contamination_process = contamination.LinearContaminationModel(T, random_state=random_state+1)
    Yt = contamination_process.sample_labels(Y)
    print("Done.")
    sys.stdout.flush()

    # Yt_batch = Yt_batch.detach().numpy()
    print("Done.")
    sys.stdout.flush()

    # Identify the set of clean observations
    conf_scores = black_box.predict_proba(X).max(axis=1)
    clean_frac = np.round(n_clean/n, decimals=5)
    threshold = np.quantile(conf_scores, 1 - clean_frac)
    I = (conf_scores >= threshold).astype(int)
    Y_obs = np.where(I == 1, Y, Yt)

    del X; torch.cuda.empty_cache()

    # Initialize an empty list to store the evaluation results
    res_list = []

    #____________________________________________________________________
    ## Estimate T using the NN algorithm
    #X_feat_torch  = torch.tensor(X_features, dtype=torch.float32)
    X_cifar_feat_torch = X_cifar_feat.to(device)
    Y_obs_torch = torch.tensor(Y_obs, dtype=torch.long)
    I_torch = torch.tensor(I, dtype=torch.long)

    if not contamination_exp_flag:
        """
        #____________________________________________________________________
        ## Estimate T using the NN with MLP
        print("Estimating T using the NN...", end=' ')
        sys.stdout.flush()
        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[64,32], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN = model_NN.contamination.contamination_matrix()
        T_hat_NN = T_hat_NN.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()
        """

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[64], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_cifar, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN cifar', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features and easier MLP...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[16,8], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_cifar, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN cifar light', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features and SLL...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_cifar, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN cifar SLL', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        """
        #____________________________________________________________________
        ## Estimate T using the NN with ResNet
        print("Estimating T using the NN with ResNet...", end=' ')
        sys.stdout.flush()
        model_resnet = NoisyLabelNet(K=K, backbone_model_="resnet", freeze_features=True, contamination_model_="uniform", epsilon_init=epsilon_init)
        history1 = train_alternate(model_resnet, X, Y_obs_torch, I_torch,
                            n_epochs=20, n_grad_steps=50,
                            batch_size=128, lr=1e-3)
        
        for param in model_resnet.backbone.net.parameters():
            param.requires_grad_(True)

        history2 = train_alternate(model_resnet, X, Y_obs_torch, I_torch,
                                n_epochs=50, n_grad_steps=50,
                               batch_size=128, lr=1e-4)

        T_hat_resnet = model_resnet.contamination.contamination_matrix()
        T_hat_resnet = T_hat_resnet.detach().numpy()

        performances = evaluate_estimate(T, T_hat_resnet, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='ResNet', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()
        """

    else:
        """
        #____________________________________________________________________
        ## Estimate T using the NN with MLP
        print("Estimating T using the NN...", end=' ')
        sys.stdout.flush()
        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[64,32], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN = model_NN.contamination.contamination_matrix()
        T_hat_NN = T_hat_NN.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()
        """

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[64], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_cifar, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN cifar gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[16,8], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_cifar, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN cifar light gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()


        #____________________________________________________________________
        ## Estimate T using the NN with cifar features and MLP
        print("Estimating T using the NN with cifar features and SLL...", end=' ')
        sys.stdout.flush()
        model_NN_cifar = NoisyLabelNet(input_dim=num_var_cifar, K=K, hidden_dims=[], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_cifar, X_cifar_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_cifar = model_NN_cifar.contamination.contamination_matrix()
        T_hat_NN_cifar = T_hat_NN_cifar.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_cifar, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN cifar SLL', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        """
        #____________________________________________________________________
        ## Estimate T using the NN algorithm with MLP and general contamination
        print("Estimating T using NN with MLP and general contamination...", end=' ')
        sys.stdout.flush()

        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[64, 32], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN = model_NN.contamination.contamination_matrix()
        T_hat_NN = T_hat_NN.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN alt gen',  **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN algorithm with ResNet and general contamination
        print("Estimating T using NN with ResNet and general contamination...", end=' ')
        sys.stdout.flush()

        model_resnet_gen = NoisyLabelNet(K=K, backbone_model_="resnet", freeze_features=True, contamination_model_="general", epsilon_init=epsilon_init)
        history1 = train_alternate(model_resnet_gen, X, Y_obs_torch, I_torch,
                            n_epochs=20, n_grad_steps=50,
                            batch_size=128, lr=1e-3)
        for param in model_resnet_gen.backbone.net.parameters():
            param.requires_grad_(True)

        history2 = train_alternate(model_resnet_gen, X, Y_obs_torch, I_torch,
                                n_epochs=50, n_grad_steps=50,
                                batch_size=128, lr=1e-4)

        T_hat_resnet_gen = model_resnet_gen.contamination.contamination_matrix()
        T_hat_resnet_gen = T_hat_resnet_gen.detach().numpy()

        performances = evaluate_estimate(T, T_hat_resnet_gen, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='ResNet gen',  **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()
        """

    # Combine all results into a single DataFrame
    res = pd.concat(res_list, ignore_index=True)
    #print(res)
    return res

# Run all experiments
results = pd.DataFrame({})
for batch in np.arange(1,num_exp+1):
    res = run_experiment(1000*seed+batch-1000)
    results = pd.concat([results, res])

    # Save results
    outfile = "results/" + outfile_prefix + ".txt"
    results.to_csv(outfile, index=False, float_format="%.5f")

print("\nFinished.\nResults written to {:s}\n".format(outfile))
sys.stdout.flush()
