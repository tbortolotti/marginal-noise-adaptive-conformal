import numpy as np
from sklearn.model_selection import train_test_split
import torch
import json
import os
import hydra
import pytorch_lightning as pl
from hydra.utils import instantiate
from hydra.core.global_hydra import GlobalHydra

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
import sys

sys.path.append("..")
sys.path.append("../third_party")

from cln import contamination
from cln.T_estimation import evaluate_estimate
from cln.T_estimation_NN import NoisyLabelNet, train_alternate
from third_party.bigearthnet.datamodules.bigearthnet_datamodule import BigEarthNetDataModule
#from third_party.bigearthnet.models.bigearthnet_module import BigEarthNetModule, BigEarthNetFeatureExtractor, TorchGeoFeatureExtractor
from third_party.bigearthnet.models.bigearthnet_module import BigEarthNetModule, TorchGeoFeatureExtractor

# Define default parameters
exp_num = 1021
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n = 1000
n_clean = 500
pi_clean = 0
contamination_exp_flag = False
seed = 1

# Parse input parameters
if True:
    print('Number of arguments:', len(sys.argv), 'arguments.')
    print('Argument List:', str(sys.argv))
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
data_name = "bigearthnet"
K = 6
num_exp = 5
epsilon_init = 0

if n_clean == 0:
    n_clean = int(np.round(pi_clean * n))
    n_noisy = n - n_clean
else:
    n_noisy = n
    n = n_noisy + n_clean

batch_size = n

device = "cuda" if torch.cuda.is_available() else "cpu"

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

## BigEarthNet setup
# Clear any previous Hydra instances
if GlobalHydra.instance().is_initialized():
    GlobalHydra.instance().clear()
hydra.initialize(config_path="../third_party/bigearthnet/configs", version_base="1.2")
cfg = hydra.compose(config_name="config")
transforms = instantiate(cfg.transforms.obj)

# Load the clean label mapping (v1 -> v2 grouped labels)
with open('../third_party/bigearthnet/data/label_mapping.json', 'r') as f:
    label_mapping = json.load(f)

dataset_name = cfg.datamodule.dataset_name
file_path_train = os.path.join(
    '../third_party/bigearthnet/data',
    f'train_{dataset_name}.csv'
)
v1v2_corresp_train = pd.read_csv(file_path_train, header=0)

# Load the pre-trained BigEarthNet model and use it as a feature extractor.
# We call it in eval mode and extract the penultimate-layer embeddings
black_box_model = BigEarthNetModule(cfg)
mod_path = os.path.join(cfg.out_directory.dir, "trained_model.pth")
black_box_model.load_state_dict(torch.load(mod_path))
black_box_model.eval()

#feature_extractor = BigEarthNetFeatureExtractor(black_box_model)
feature_extractor = TorchGeoFeatureExtractor(black_box_model)

# Define dataloader
datamodule = BigEarthNetDataModule(cfg.datamodule.dataset_dir,
                                   cfg.datamodule.dataset_name,
                                   batch_size,
                                   cfg.datamodule.num_workers,
                                   transforms,
                                   label_mapping)
datamodule.setup()

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


def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    print("\nGenerating data...", end=' ')
    sys.stdout.flush()
    dataloader = datamodule.train_dataloader(seed=random_state)
    batch = next(iter(dataloader))
    X = batch['data']

    print("Shape:", X.shape)
    print("Min:", X.min().item(), "Max:", X.max().item())

    # Retrieve the corresponding clean (v2-grouped) labels from the CSV,
    # using the same generator the datamodule used to shuffle its dataset.
    shuffled_csv = v1v2_corresp_train.iloc[datamodule.last_train_indices].reset_index(drop=True)
    batch_csv = shuffled_csv.iloc[:batch_size]
    Y = batch_csv['v2-labels-grouped'].to_numpy()

    # Drop rows where clean label is NaN
    valid = ~np.isnan(Y)
    X = X[valid]
    Y = Y[valid].astype(int)
    print(f"Done. Batch size after NaN removal: {len(Y)}")
    sys.stdout.flush()

    # Feature extraction
    X_feat = feature_extractor.transform(X).numpy()
    num_var = X_feat.shape[1]

    # Artificial label contamination
    print("Generating contaminated labels...", end=' ')
    sys.stdout.flush()
    contamination_process = contamination.LinearContaminationModel(T, random_state=random_state + 1)
    Yt = contamination_process.sample_labels(Y)
    print("Done.")
    sys.stdout.flush()

    # Identify the set of clean observations
    print("Generating clean dataset...", end=' ')
    sys.stdout.flush()
    # Identify the central observations in the training set to build the clean dataset
    conf_scores = black_box_model.predict_proba(X).max(axis=1)
    top_indices = np.argsort(conf_scores)[-n_clean:]
    I = np.zeros(len(Y))
    I[top_indices] = 1
    Y_obs = np.where(I == 1, Y, Yt)
    print("Done.")
    sys.stdout.flush()
    del X; torch.cuda.empty_cache()

    # Initialize an empty list to store the evaluation results
    res_list = []

    #____________________________________________________________________
    ## Estimate T using the NN algorithm
    #X_feat_torch  = torch.tensor(X_features, dtype=torch.float32)
    X_feat_torch = torch.tensor(X_feat, dtype=torch.float32).to(device)
    Y_obs_torch = torch.tensor(Y_obs, dtype=torch.long)
    I_torch = torch.tensor(I, dtype=torch.long)

    if not contamination_exp_flag:
        #____________________________________________________________________
        ## Estimate T using the NN with features and MLP
        print("Estimating T using the NN with features...", end=' ')
        sys.stdout.flush()
        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[64], contamination_model_="uniform", epsilon_init=epsilon_init)
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

        #____________________________________________________________________
        ## Estimate T using the NN with features and MLP
        print("Estimating T using the NN with features and easier MLP...", end=' ')
        sys.stdout.flush()
        model_NN_light = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN_light, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_light, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_light = model_NN_light.contamination.contamination_matrix()
        T_hat_NN_light = T_hat_NN_light.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_light, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN light', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with features and MLP
        print("Estimating T using the NN with features and SLL...", end=' ')
        sys.stdout.flush()
        model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
        train_alternate(model_NN_sll, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_sll, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_sll = model_NN_sll.contamination.contamination_matrix()
        T_hat_NN_sll = T_hat_NN_sll.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_sll, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN SLL', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

    else:
        #____________________________________________________________________
        ## Estimate T using the NN with features and MLP
        print("Estimating T using the NN with features...", end=' ')
        sys.stdout.flush()
        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[64], contamination_model_="general", epsilon_init=epsilon_init)
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

        #____________________________________________________________________
        ## Estimate T using the NN with features and MLP
        print("Estimating T using the NN with features...", end=' ')
        sys.stdout.flush()
        model_NN_light = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_light, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_light, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_light = model_NN_light.contamination.contamination_matrix()
        T_hat_NN_light = T_hat_NN_light.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_light, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN light gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with features and MLP
        print("Estimating T using the NN with features and SLL...", end=' ')
        sys.stdout.flush()
        model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_sll, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-2, verbose=False)
        train_alternate(model_NN_sll, X_feat_torch, Y_obs_torch, I_torch, n_epochs=50, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_sll= model_NN_sll.contamination.contamination_matrix()
        T_hat_NN_sll = T_hat_NN_sll.detach().numpy()

        performances = evaluate_estimate(T, T_hat_NN_sll, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN SLL gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

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