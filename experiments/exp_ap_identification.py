import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.dummy import DummyClassifier

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm
import pdb
import copy

import sys
sys.path.append("..")
sys.path.append("../third_party")


from cln import data
from cln import contamination
from cln.AP_identification import AnchorPointsIdentification
from cln.T_estimation import evaluate_estimate, TMatrixEstimation
from third_party import arc

# Define default parameters
exp_num = 801

# Data simulation parameters
data_name = 'syntheticAP'
num_var = 2
K = 4
pi_easy=1
delta_shift=0
center_scale=1
sigma_easy = 1
sigma_hard = 1
R_easy = 1.5
R_hard = 4.0
flipy = 0

# Point predictor
model_name = 'RFC'

# Contamination parameters
epsilon = 0.2
contamination_model = "uniform"

# Set sizes
n_train1 = 10000
n_train2 = 5000
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 19:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()
    exp_num = int(sys.argv[1])
    data_name = sys.argv[2]
    num_var = int(sys.argv[3])
    K = int(sys.argv[4])
    pi_easy = float(sys.argv[5])
    delta_shift = float(sys.argv[6])
    center_scale = float(sys.argv[7])
    sigma_easy = float(sys.argv[8])
    sigma_hard = float(sys.argv[9])
    R_easy = float(sys.argv[10])
    R_hard = float(sys.argv[11])
    flipy = float(sys.argv[12])
    model_name = sys.argv[13]
    epsilon = float(sys.argv[14])
    contamination_model = sys.argv[15]
    n_train1 = int(sys.argv[16])
    n_train2 = int(sys.argv[17])
    seed = int(sys.argv[18])

# Define other parameters
A_easy = sigma_easy * np.eye(2)
A_hard = sigma_hard * np.eye(2)
nu = 0

batch_size = 10
gamma_vec = np.asarray([0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2], dtype=float)

# Initialize the data distribution
if data_name == "syntheticAP":
    data_distribution = data.DataModel_AP(K, num_var,
                                          pi_easy=pi_easy, delta_shift=delta_shift, center_scale=center_scale, epsilon0=flipy,
                                          distribution_type="truncated_gaussian",
                                          A_easy=A_easy, R_easy=R_easy, A_hard=A_hard, R_hard=R_hard,
                                          random_state=seed)
else:
    print("Unknown data distribution!")
    sys.stdout.flush()
    exit(-1)

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
if model_name == 'RFC':
    black_box = arc.black_boxes.RFC(n_estimators=100, max_features="sqrt")
elif model_name == 'SVC':
    black_box = arc.black_boxes.SVC(clip_proba_factor = 1e-5)
elif model_name == 'NN':
    black_box = arc.black_boxes.NNet(max_iter=100)
elif model_name == 'MF':
    black_box = DummyClassifier(strategy="most_frequent")
elif model_name == 'SRC':
    black_box = DummyClassifier(strategy="stratified")
else:
    print("Unknown model!")
    sys.stdout.flush()
    exit(-1)

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'num_var':[num_var], 'K':[K],
                       'pi_easy':[pi_easy], 'delta_shift':[delta_shift], 'center_scale':[center_scale],
                       'sigma_easy':[sigma_easy], 'sigma_hard':[sigma_hard], 'R_easy':[R_easy], 'R_hard':[R_hard],
                       'n_train1':[n_train1], 'n_train2':[n_train2],
                       'flipy':[flipy], 'epsilon':[epsilon], 'contamination':[contamination_model],
                       'model_name':[model_name],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_p" + str(num_var) + "_K" + str(K)
outfile_prefix += "_pi_easy" + str(pi_easy) + "_dshift" + str(delta_shift) + "_cscale" + str(center_scale)
outfile_prefix += "_seasy" + str(sigma_easy) + "_shard" + str(sigma_hard) + "_Reasy" + str(R_easy) + "_Rhard" + str(R_hard)
outfile_prefix += "_flipy" + str(flipy)
outfile_prefix += "_"  + model_name
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
    data_distribution.set_seed(random_state+1)
    X, Y = data_distribution.sample(n_train1+n_train2)
    print("Done.")
    sys.stdout.flush()

    # Generate the contaminated labels
    print("Generating contaminated labels...", end=' ')
    sys.stdout.flush()
    contamination_process = contamination.LinearContaminationModel(T, random_state=random_state+2)
    Yt = contamination_process.sample_labels(Y)
    print("Done.")
    sys.stdout.flush()

    # Separate data into training and calibration
    X_train1, X_train2, Y_train1, Y_train2, Yt_train1, Yt_train2 = train_test_split(X, Y, Yt, test_size=n_train2, random_state=random_state+3)

    methods = {
        "benchmark": lambda: AnchorPointsIdentification(X_train1, Yt_train1, X_train2, Yt_train2, K,
                                            use_classifier=True, black_box=black_box,
                                            calibrate_gamma=True),

        "EE": lambda: AnchorPointsIdentification(X_train1, Yt_train1, X_train2, Yt_train2, K,
                                            outlier_detection=True,
                                            outlier_detection_method="elliptic_envelope",
                                            selection="accuracy"),

        "IF": lambda: AnchorPointsIdentification(X_train1, Yt_train1, X_train2, Yt_train2, K,
                                                    outlier_detection=True,
                                                    outlier_detection_method="isolation_forest",
                                                    selection="accuracy"),

        "LOF": lambda: AnchorPointsIdentification(X_train1, Yt_train1, X_train2, Yt_train2, K,
                                                    outlier_detection=True,
                                                    outlier_detection_method="lof",
                                                    selection="accuracy")
    }

    # Initialize an empty list to store the evaluation results
    res = pd.DataFrame({})
    res_list = []

    # Estimate using all the clean/noisy labels correspondence
    T_method = TMatrixEstimation(X_train2, Y_train2, Yt_train2, K, estimation_method="empirical")
    T_hat = T_method.get_estimate()
    
    performances = evaluate_estimate(T, T_hat, Y_train2, Y_train2, Yt_train2, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='Clean sample',  gamma_opt=1, **performances)
    res_list.append(res_update)

    # Loop through the methods, apply them, and evaluate the results
    for method_name, method_func in methods.items():
        print(f"Applying {method_name} method...", end=' ')
        sys.stdout.flush()

        # Initialize and apply the method for anchor points
        method = method_func()
        Ya_train2, _, _, _ = method.get_anchor_points()

        # Use anchor points to estimate T
        T_method = TMatrixEstimation(X_train2, Ya_train2, Yt_train2, K, estimation_method="empirical_parametricRR")
        T_hat = T_method.get_estimate()

        print("Done.")
        sys.stdout.flush()

        performances = evaluate_estimate(T, T_hat, Y_train2, Ya_train2, Yt_train2, K, epsilon0=flipy)
        res_update = header.copy()
        res_update = res_update.assign(Method=method_name, n_train1=n_train1, n_train2=n_train2, **performances)
        res_list.append(res_update)

    # Combine all results into a single DataFrame
    res = pd.concat(res_list, ignore_index=True)
    #print(res)
    return res

# Run all experiments
results = pd.DataFrame({})
for batch in np.arange(1,batch_size+1):
    res = run_experiment(1000*seed+batch-1000)
    results = pd.concat([results, res])

    # Save results
    outfile = "results/" + outfile_prefix + ".txt"
    results.to_csv(outfile, index=False, float_format="%.5f")

print("\nFinished.\nResults written to {:s}\n".format(outfile))
sys.stdout.flush()
