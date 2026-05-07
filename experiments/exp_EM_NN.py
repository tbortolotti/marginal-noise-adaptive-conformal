import numpy as np
from sklearn.model_selection import train_test_split
import torch

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm
import pdb
import copy
from sklearn.ensemble import RandomForestClassifier

import sys
sys.path.append("..")
sys.path.append("../third_party")


from cln import data
from cln import contamination
from cln.T_estimation import TMatrixEstimation
from cln.T_estimation_EM import Dataset, run_em, predict
from cln.T_estimation_NN import NoisyLabelNet, train, train_alternate, train_em_style
from cln.utils import evaluate_predictions, estimate_rho
from cln.classification import MarginalLabelNoiseConformal
from third_party import arc


# Define default parameters
exp_num = 711
data_name = 'synthetic6'
model_name = 'RFC'
num_var = 20
K = 4
epsilon = 0.2
contamination_model = "uniform"
n_train = 10000
n_clean = 500
pi_clean = 0
n_cal = 5000
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 13:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()
    exp_num = int(sys.argv[1])
    model_name = sys.argv[2]
    data_name = sys.argv[3]
    num_var = int(sys.argv[4])
    K = int(sys.argv[5])
    epsilon = float(sys.argv[6])
    contamination_model = sys.argv[7]
    n_train = int(sys.argv[8])
    n_clean = int(sys.argv[9])
    pi_clean = float(sys.argv[10])
    n_cal = int(sys.argv[11])
    seed = int(sys.argv[12])

# Define other constant parameters
estimate = "none"
n_test = 2000
batch_size = 10
allow_empty = True
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000
nu = 0.2

epsilon_init = 0

if n_clean == 0:
    n_clean = int(np.round(pi_clean * n_train))

# Initialize the data distribution
if data_name == "synthetic1":
    data_distribution = data.DataModel_1(K, num_var, signal=1, random_state=seed)
elif data_name == "synthetic2":
    data_distribution = data.DataModel_2(K, num_var, signal=1, random_state=seed)
elif data_name == "synthetic3":
    data_distribution = data.DataModel_3(K, num_var, signal=1, random_state=seed)
elif data_name == "synthetic6":
    data_distribution = data.DataModel_6(K, num_var, random_state=seed)
else:
    print("Unknown data distribution!")
    sys.stdout.flush()
    exit(-1)

# Estimate the label proportions from the population model
rho = data_distribution.estimate_rho()

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
elif contamination_model=="mild":
    # --- Option 1: mild, nearly diagonal contamination ---
    # Each label is mostly correct, small symmetric noise
    T = np.array([
        [0.80, 0.10, 0.05, 0.05],
        [0.10, 0.80, 0.05, 0.05],
        [0.05, 0.05, 0.80, 0.10],
        [0.05, 0.05, 0.10, 0.80],
    ], dtype=np.float32)
elif contamination_model=="asymmetric":
    # --- Option 2: asymmetric contamination ---
    # Label 0 often gets flipped to label 1, but not vice versa
    # Labels 2 and 3 are relatively clean
    T = np.array([
        [0.60, 0.05, 0.05, 0.05],
        [0.30, 0.75, 0.05, 0.10],
        [0.05, 0.10, 0.80, 0.10],
        [0.05, 0.10, 0.10, 0.75],
    ], dtype=np.float32)
elif contamination_model=="hard":
    # --- Option 3: hard, nearly uniform contamination ---
    # Very little signal left in the contaminated labels
    T = np.array([
        [0.40, 0.20, 0.20, 0.20],
        [0.20, 0.40, 0.20, 0.20],
        [0.20, 0.20, 0.40, 0.20],
        [0.20, 0.20, 0.20, 0.40],
    ], dtype=np.float32)
else:
    print("Unknown contamination (M) model!")
    sys.stdout.flush()
    exit(-1)

# Compute the contaminated label proportions
rho_tilde = np.dot(T, rho)

# Initialize black-box model
if model_name == 'RFC':
    black_box = arc.black_boxes.RFC(n_estimators=100)
elif model_name == 'SVC':
    black_box = arc.black_boxes.SVC(clip_proba_factor = 1e-5)
elif model_name == 'NN':
    black_box = arc.black_boxes.NNet(max_iter=100)
else:
    print("Unknown model!")
    sys.stdout.flush()
    exit(-1)

# Initialize black-box model
black_box_SVC = arc.black_boxes.SVC(clip_proba_factor = 1e-5)

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'num_var':[num_var], 'K':[K],
                       'n_train':[n_train], 'n_clean':[n_clean], 'n_cal':[n_cal],
                       'epsilon':[epsilon], 'contamination':[contamination_model],
                       'model_name':[model_name], 'estimate':[estimate], 'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_p" + str(num_var)
outfile_prefix += "_K" + str(K) + "_" + model_name
outfile_prefix += "_eps" + str(epsilon) + "_" + contamination_model
outfile_prefix += "_nt_" + str(n_train) + "_ncl_" + str(n_clean) + "_picl_" + str(pi_clean) +"_nc" + str(n_cal) + "_est" + estimate + "_seed" + str(seed)
print("Output file: {:s}.".format("results/"+outfile_prefix), end="\n")
sys.stdout.flush()

# Describe the experiment
def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    # Generate a large data set
    print("\nGenerating data...", end=' ')
    sys.stdout.flush()
    data_distribution.set_seed(random_state+1)
    X_all, Y_all = data_distribution.sample(n_train+n_cal+n_test)
    print("Done.")
    sys.stdout.flush()

    # Separate the test set
    X, X_test, Y, Y_test = train_test_split(X_all, Y_all, test_size=n_test, random_state=random_state+2)

    # Generate the contaminated labels
    print("Generating contaminated labels...", end=' ')
    sys.stdout.flush()
    contamination_process = contamination.LinearContaminationModel(T, random_state=random_state+3)
    Yt = contamination_process.sample_labels(Y)
    print("Done.")
    sys.stdout.flush()

    # Estimate (if applicable) the label contamination model
    if estimate=="none":
        rho_tilde_hat = rho_tilde
    elif estimate=="rho":
        rho_tilde_hat = estimate_rho(Yt, K)
    else:
        print("Unknown estimation option!")
        sys.stdout.flush()
        exit(-1)

    # Separate data into training and calibration
    X_train, X_cal, Y_train, Y_cal, Yt_train, Yt_cal = train_test_split(X, Y, Yt, test_size=n_cal, random_state=random_state+4)

    #_________________________________________________________________
    # Generate the clean dataset
    # This is done by relying on two independent datasets, one where I build
    # a classifier and the other where I use the classifier and select the
    # n_clean observations that are the least difficult to predict
    n_train0 = 10000
    data_distribution.set_seed(random_state+5)
    X_train0, Y_train0 = data_distribution.sample(n_train0)
    Yt_train0 = contamination_process.sample_labels(Y_train0)
    # Fit a quick RFC on noisy labels to score each observation's difficulty
    rfc_easy = RandomForestClassifier(n_estimators=100, random_state=random_state+4)
    rfc_easy.fit(X_train0, Yt_train0)

    # Assign I=1 to the top clean_frac fraction by confidence
    conf_scores = rfc_easy.predict_proba(X_train).max(axis=1)
    clean_frac = np.round(n_clean/n_train, decimals=5)
    threshold = np.quantile(conf_scores, 1 - clean_frac)
    I_train = (conf_scores >= threshold).astype(int)

    X_clean = X_train[I_train==1]
    Y_clean = Y_train[I_train==1]
    Yt_clean = Yt_train[I_train==1]
    Yt_train[I_train==1] = Y_clean
    #_________________________________________________________________

    # 
    # Fit the point predictor on the training set (with also the clean samples)
    black_box_pt = copy.deepcopy(black_box)
    black_box_pt.fit(X_train, Yt_train)

    # Estimate T using all the clean/noisy labels correspondence on clean dataset
    T_method = TMatrixEstimation(Y_clean, Yt_clean, K, estimation_method="empirical_parametricRR")
    T_hat_clean = T_method.get_estimate()

    # EM method for T estimation
    print("Estimating T using EM algorithm...", end=' ')
    sys.stdout.flush()
    X_intercept = np.hstack([np.ones((n_train, 1)), X_train])
    data = Dataset(X=X_intercept, Y_obs=Yt_train, I=I_train, K=K)
    result_EM = run_em(data, contamination_model_="uniform", eps_init=epsilon_init, max_iter=100, tol=1e-7, verbose=False)
    T_hat_EM = result_EM.T
    print("Done.")
    sys.stdout.flush()


    #____________________________________________________________________
    X_torch  = torch.tensor(X_train, dtype=torch.float32)
    Yt_torch = torch.tensor(Yt_train, dtype=torch.long)
    I_torch = torch.tensor(I_train, dtype=torch.long)

    #____________________________________________________________________
    ## Estimate T using the NN algorithm with SLL
    print("Estimating T using the NN with SLL...", end=' ')
    sys.stdout.flush()
    model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
    train_alternate(model_NN_sll, X_torch, Yt_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=5e-2, verbose=False)
    train_alternate(model_NN_sll, X_torch, Yt_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

    T_hat_NN_sll = model_NN_sll.contamination.contamination_matrix()
    T_hat_NN_sll = T_hat_NN_sll.detach().numpy()
    print("Done.")
    sys.stdout.flush()

     #____________________________________________________________________
    ## Estimate T using the NN and alternate training
    print("Estimating T using the NN...", end=' ')
    sys.stdout.flush()
    model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="uniform", epsilon_init=epsilon_init)
    train_alternate(model_NN, X_torch, Yt_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=5e-2, verbose=False)
    train_alternate(model_NN, X_torch, Yt_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

    T_hat_NN = model_NN.contamination.contamination_matrix()
    T_hat_NN = T_hat_NN.detach().numpy()
    print("Done.")
    sys.stdout.flush()

    alpha = 0.1
    guarantee = 'marginal'

    res = pd.DataFrame({})

    print("\nSeeking {:s} coverage at level {:.2f}.".format(guarantee, 1-alpha))

    # Define a dictionary of methods with their names and corresponding initialization parameters
    methods = {
        "Standard": lambda: arc.methods.SplitConformal(X_cal, Yt_cal, black_box_pt, K, alpha, n_cal=-1,
                                                       pre_trained=True, random_state=random_state),

        "Standard using clean": lambda: arc.methods.SplitConformal(X_clean, Y_clean, black_box_pt, K, alpha, n_cal=-1,
                                                                pre_trained=True, random_state=random_state),

        "Adaptive optimized+": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box_pt, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ clean": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box_pt, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_clean, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ EM": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box_pt, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_EM, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ NN SLL": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box_pt, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_NN_sll, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

        "Adaptive optimized+ NN": lambda: MarginalLabelNoiseConformal(X_cal, Yt_cal, black_box_pt, K, alpha, n_cal=-1,
                                                                    epsilon=epsilon, T=T_hat_NN, rho_tilde=rho_tilde_hat,
                                                                    allow_empty=allow_empty, method="improved",
                                                                    optimized=True, optimistic=True, verbose=False,
                                                                    pre_trained=True, random_state=random_state),

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
for batch in np.arange(1,batch_size+1):
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
