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
from cln.AP_identification import AnchorPointsIdentification
from cln.T_estimation import TMatrixEstimation
from cln.utils import evaluate_predictions, estimate_rho
from cln.classification import MarginalLabelNoiseConformal
from third_party import arc

from data_torch import Cifar10DataSet, ImageNetResNet18Features, ResNet18
from torch.utils.data import DataLoader

import gc
import copy


# Define default parameters
exp_num = 1001
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n_train1 = 1000
n_train2 = 1000
n_cal = 1000
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
    n_cal = int(sys.argv[7])
    seed = int(sys.argv[8])

# Define other constant parameters
data_name = "cifar10"
K = 10
num_exp = 5
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000
n_test = 500
batch_size = n_train1 + n_train2 + n_cal + n_test

allow_empty = True
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000

# Set default directories
data_dir = "/home1/tb_214/data/cifar10"
noisy_data_dir = "/home1/tb_214/data/cifar10h"

print(f"Data Directory: {data_dir}")
print(f"Noisy Data Directory: {noisy_data_dir}")

dataset = Cifar10DataSet(data_dir=data_dir, noisy_data_dir=noisy_data_dir, random_state=2024)

#print(f"\nOverall numerosity of dataset {len(dataset)}")
#sys.stdout.flush()

loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, num_workers=1)
feature_extractor = ImageNetResNet18Features()

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
    print("Unknown contamination (M) model!")
    sys.stdout.flush()
    exit(-1)

# Initialize black-box model
black_box = ResNet18()
#black_box = arc.black_boxes.MLPBlackBox(clip_proba_factor=1e-5)

# Initialize the black-box model
black_box_SVC = arc.black_boxes.SVC(clip_proba_factor = 1e-5)

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'K':[K],
                       'n_train1':[n_train1], 'n_train2':[n_train2], 'n_cal':[n_cal],
                       'epsilon':[epsilon], 'contamination':[contamination_model],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_n" + str(batch_size)
outfile_prefix += "_eps" + str(epsilon) + "_" + contamination_model
outfile_prefix += "_nt1_" + str(n_train1) + "_nt2_" + str(n_train2) +"_nc" + str(n_cal) + "_seed" + str(seed)
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
    X_all, X_all_imagenet, Y_all, _, _, _, _ = next(iter(loader))
    Y_all = Y_all.detach().numpy()
    print("Done.")
    sys.stdout.flush()

    # Separate the test set
    X, X_test, X_imagenet, _, Y, Y_test = train_test_split(X_all, X_all_imagenet, Y_all, test_size=n_test, random_state=random_state+1)
    del X_all, X_all_imagenet, Y_all

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
    print(rho_tilde_hat)
    print("Done.")
    sys.stdout.flush()


    # Separate data into training and calibration
    X_imagenet_train, X_imagenet_cal, _, X_cal, Y_train, Y_cal, Yt_train, Yt_cal = train_test_split(X_imagenet, X, Y, Yt, test_size=n_cal, random_state=random_state+3)
    del X_imagenet, X, Y, Yt

    # Extract features
    print("Extract features for T estimation...", end=' ')
    sys.stdout.flush()
    # Split training set in two, as it'll be needed to estimate T
    X_train1, X_train2, _, Y_train2, Yt_train1, Yt_train2 = train_test_split(X_imagenet_train, Y_train, Yt_train, test_size=n_train2, random_state=random_state+4)
    del X_imagenet_train, Y_train, Yt_train

    # Operate transformation of X to fit SVC and identify anchor points
    X_features_train1 = feature_extractor.transform(X_train1).numpy()
    del X_train1; torch.cuda.empty_cache()

    X_features_train2 = feature_extractor.transform(X_train2).numpy()
    del X_train2; torch.cuda.empty_cache()

    X_features_cal = feature_extractor.transform(X_imagenet_cal).numpy()
    del X_imagenet_cal; torch.cuda.empty_cache()

    print("Done.")
    sys.stdout.flush()

    # Train MLP black box on the training set
    #print("Training MLP black box...", end=' ')
    #sys.stdout.flush()
    #X_features_train_all = np.concatenate([X_features_train1, X_features_train2])
    #Y_train_all = np.concatenate([Yt_train1, Yt_train2])

    #black_box_MLP = copy.deepcopy(black_box)
    #black_box_MLP.fit(X_features_train_all, Y_train_all)
    #del X_features_train_all, Y_train_all
    #torch.cuda.empty_cache()
    #print("Done.")
    #sys.stdout.flush()

    print("Estimating contamination matrix...", end=' ')
    sys.stdout.flush()
    # Estimate T using all the clean/noisy labels correspondence
    T_method = TMatrixEstimation(Y_train2, Yt_train2, K, estimation_method="empirical_parametricRR")
    T_hat_clean = T_method.get_estimate()

    ## Anchor points method for T estimation, using SVC as classifier
    method = AnchorPointsIdentification(X_features_train1, Yt_train1, X_features_train2, Yt_train2, K,
                                        use_classifier=True, black_box=black_box_SVC,
                                        calibrate_gamma=True)
    Ya_train2, _, _, _ = method.get_anchor_points()
    T_method = TMatrixEstimation(Ya_train2, Yt_train2, K, estimation_method="empirical_parametricRR")
    T_hat_SVC = T_method.get_estimate()

    # Anchor points method for T estimation, using combination of outlier detectors as classifier
    method = AnchorPointsIdentification(X_features_train1, Yt_train1, X_features_train2, Yt_train2, K,
                                        black_box=black_box_SVC,
                                        optimal_method=True,
                                        random_state=random_state+5)
    Ya_train2, _, _, _ = method.get_anchor_points()
    T_method = TMatrixEstimation(Ya_train2, Yt_train2, K, estimation_method="empirical_parametricRR")
    T_hat_opt = T_method.get_estimate()

    del X_features_train2
    print("Done.")
    sys.stdout.flush()

    print("Identifying set of anchor points...", end=' ')
    sys.stdout.flush()
    # Create dataset of sole anchor points
    method = AnchorPointsIdentification(X_features_train1, Yt_train1, X_features_cal, Yt_cal, K,
                                        black_box=black_box_SVC,
                                        optimal_method=True,
                                        random_state=random_state+6)
    Ya_cal, _, _, _ = method.get_anchor_points()
    idxs_anchor = (Ya_cal != -1)
    X_anchor = X_cal[idxs_anchor,]
    Y_anchor = Ya_cal[idxs_anchor]
    del X_features_train1, X_features_cal
    #Y_anchor = Y_cal[idxs_anchor]
    print("Done.")
    sys.stdout.flush()

    # Compute size of AP set and imbalance ratio
    size_ap = np.sum(Y_anchor!=-1)
    Y_anchor_valid = Y_anchor[Y_anchor != -1]
    counts = np.bincount(Y_anchor_valid, minlength=K)
    IR = np.inf if np.any(counts == 0) else counts.max() / counts.min()

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

        "Standard using AP": lambda: arc.methods.SplitConformal(X_anchor, Y_anchor, black_box, K, alpha, n_cal=-1,
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

        "Adaptive optimized+ AP SVC": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_SVC, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ AP opt": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_opt, rho_tilde=rho_tilde_hat,
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
        res_new['size_ap'] = size_ap
        res_new['ir'] = IR
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
