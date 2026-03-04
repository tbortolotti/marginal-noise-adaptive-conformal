import numpy as np
from sklearn.model_selection import train_test_split
import torch
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

from cln import contamination
from cln.AP_identification import AnchorPointsIdentification
from cln.T_estimation import evaluate_estimate, TMatrixEstimation
from third_party import arc

from data_torch import Cifar10DataSet, ResNet18
from torch.utils.data import DataLoader

# Define default parameters
exp_num = 901
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n_train1 = 1000
n_train2 = 1000
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 9:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()
    exp_num = int(sys.argv[1])
    epsilon = float(sys.argv[2])
    nu = float(sys.argv[3])
    contamination_model = sys.argv[4]
    n_train1 = int(sys.argv[5])
    n_train2 = int(sys.argv[6])
    seed = int(sys.argv[7])


# Define other constant parameters
batch_size = n_train1 + n_train2
data_name = "cifar10"
K = 10
num_exp = 1
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000

# Set default directories
data_dir = "/home/tb_214/data/cifar10"
noisy_data_dir = "/home/tb_214/data/cifar10h"

print(f"Data Directory: {data_dir}")
print(f"Noisy Data Directory: {noisy_data_dir}")

dataset = Cifar10DataSet(data_dir, noisy_data_dir, normalize=True, random_state=2026)
#dataset = Cifar10DataSet(data_dir, normalize=True, random_state=2026)
loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=1)

# Initialize the black-box model
black_box_RN18 = ResNet18()
black_box_SVC = arc.black_boxes.SVC(clip_proba_factor = 1e-5)

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
                       'n_train1':[n_train1], 'n_train2':[n_train2],
                       'epsilon':[epsilon], 'contamination':[contamination_model],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_n" + str(batch_size)
outfile_prefix += "_eps" + str(epsilon) + "_" + contamination_model
outfile_prefix += "_nt1_" + str(n_train1) + "_nt2_" + str(n_train2) + "_seed" + str(seed)
print("Output file: {:s}.".format("results/"+outfile_prefix), end="\n")
sys.stdout.flush()

# Describe the experiment
def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    # Generate a large data set
    print("\nGenerating data...", end=' ')
    sys.stdout.flush()
    #X_batch, Y_batch, Y_lab_batch, Yt_batch, Yt_lab_batch, idx_batch = next(iter(loader))
    X, Y, Y_lab, _, _, idx_batch = next(iter(loader))
    
    X_detach = X.detach().cpu().numpy().reshape(X.shape[0], -1)
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

    # Separate data into first-stage training and anchor-selection set
    X_train1, X_train2, X_train1_detach, X_train2_detach, Y_train1, Y_train2, Yt_train1, Yt_train2 = train_test_split(X, X_detach, Y, Yt, test_size=n_train2, random_state=random_state+2)

    methods = {
        "SVC": lambda: AnchorPointsIdentification(X_train1_detach, Yt_train1, X_train2_detach, Yt_train2, K,
                                                    use_classifier=True, black_box=black_box_SVC,
                                                    calibrate_gamma=True),

        "ResNet18": lambda: AnchorPointsIdentification(X_train1, Yt_train1, X_train2, Yt_train2, K,
                                                    use_classifier=True, black_box=black_box_RN18,
                                                    pretrained=True,
                                                    calibrate_gamma=True),

        "EE": lambda: AnchorPointsIdentification(X_train1_detach, Yt_train1, X_train2_detach, Yt_train2, K,
                                            outlier_detection=True,
                                            outlier_detection_method="elliptic_envelope",
                                            selection="accuracy"),

        "IF": lambda: AnchorPointsIdentification(X_train1_detach, Yt_train1, X_train2_detach, Yt_train2, K,
                                                    outlier_detection=True,
                                                    outlier_detection_method="isolation_forest",
                                                    selection="accuracy")
    }

    # Initialize an empty list to store the evaluation results
    res = pd.DataFrame({})
    res_list = []

    # Clean sample - Estimate using all the clean/noisy labels correspondence
    T_method = TMatrixEstimation(Y_train2, Yt_train2, K, estimation_method="empirical")
    T_hat = T_method.get_estimate()
    
    performances = evaluate_estimate(T, T_hat, Y_train2, Y_train2, Yt_train2, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='Clean sample', n_train1=n_train1, n_train2=n_train2, size=np.sum(Y_train2!=-1), **performances)
    res_list.append(res_update)

    # Loop through the methods, apply them, and evaluate the results
    for method_name, method_func in methods.items():
        print(f"Applying {method_name} method...", end=' ')
        sys.stdout.flush()

        # Initialize and apply the method for anchor points
        method = method_func()
        Ya_train2, _, _, _ = method.get_anchor_points()

        # Use anchor points to estimate T
        T_method = TMatrixEstimation(Ya_train2, Yt_train2, K, estimation_method="empirical_parametricRR")
        T_hat = T_method.get_estimate()

        print("Done.")
        sys.stdout.flush()

        performances = evaluate_estimate(T, T_hat, Y_train2, Ya_train2, Yt_train2, K)
        res_update = header.copy()
        res_update = res_update.assign(Method=method_name, n_train1=n_train1, n_train2=n_train2, size=np.sum(Ya_train2!=-1), **performances)
        res_list.append(res_update)

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
