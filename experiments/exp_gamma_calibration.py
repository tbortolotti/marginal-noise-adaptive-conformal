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
data_name = 'synthetic1'
num_var = 20
K = 4
signal = 1
model_name = 'RFC'
flipy = 0.01
epsilon = 0.2
contamination_model = "uniform"
n_train = 10000
n = 5000
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
    data_name = sys.argv[2]
    num_var = int(sys.argv[3])
    K = int(sys.argv[4])
    signal = float(sys.argv[5])
    model_name = sys.argv[6]
    flipy = float(sys.argv[7])
    epsilon = float(sys.argv[8])
    contamination_model = sys.argv[9]
    n_train = int(sys.argv[10])
    n = int(sys.argv[11])
    seed = int(sys.argv[12])

# Define other constant parameters
nu = 0
batch_size = 20
gamma_vec = np.asarray([0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2], dtype=float)

# Initialize the data distribution
if data_name == "synthetic1":
    data_distribution = data.DataModel_1(K, num_var, signal=signal, random_state=seed)
elif data_name == "synthetic1_easy":
    data_distribution = data.DataModel_1_easy(K, num_var, signal=signal, flipy=flipy, random_state=seed)
elif data_name == "synthetic2":
    data_distribution = data.DataModel_2(K, num_var, signal=signal, random_state=seed)
elif data_name == "synthetic3":
    data_distribution = data.DataModel_3(K, num_var, signal=signal, random_state=seed)
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
                       'signal':[signal], 'n_train':[n_train], 'n':[n],
                       'flipy':[flipy], 'epsilon':[epsilon], 'contamination':[contamination_model],
                       'model_name':[model_name],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_p" + str(num_var)
outfile_prefix += "_K" + str(K) + "_signal" + str(signal) + "_" + model_name
outfile_prefix += "_flipy" + str(flipy) + "_eps" + str(epsilon) + "_" + contamination_model
outfile_prefix += "_nt" + str(n_train) + "_n" + str(n) + "_seed" + str(seed)
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
    X, Y = data_distribution.sample(n_train+n)
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
    X_train1, X_train2, Y_train1, Y_train2, Yt_train1, Yt_train2 = train_test_split(X, Y, Yt, test_size=n, random_state=random_state+3)

    # Fit the point predictor on the training set
    black_box_pt = copy.deepcopy(black_box)
    black_box_pt.fit(X_train1, Yt_train1)

    methods = {
        "D2L": lambda: AnchorPointsIdentification(X_train2, Yt_train2, K, black_box_pt, calibrate_gamma=True, gamma_vec=gamma_vec, elbow_detection_method="D2L"),

        "threshold": lambda: AnchorPointsIdentification(X_train2, Yt_train2, K, black_box_pt, gamma=(50*K/n)),

        "top3perc": lambda: AnchorPointsIdentification(X_train2, Yt_train2, K, black_box_pt, gamma=0.03),

        "D2L filtered": lambda: AnchorPointsIdentification(X_train2, Yt_train2, K, black_box_pt, calibrate_gamma=True, gamma_vec=gamma_vec, elbow_detection_method="D2L", ap_filter=True),

        "threshold filtered": lambda: AnchorPointsIdentification(X_train2, Yt_train2, K, black_box_pt, gamma=(50*K/n), ap_filter=True),

        "top3perc filtered": lambda: AnchorPointsIdentification(X_train2, Yt_train2, K, black_box_pt, gamma=0.03, ap_filter=True),
    }

    # Initialize an empty list to store the evaluation results
    res = pd.DataFrame({})
    res_list = []

    # Estimate using all the clean/noisy labels correspondence
    T_method = TMatrixEstimation(X_train2, Y_train2, Yt_train2, K, estimation_method="empirical")
    T_hat = T_method.get_estimate()
    
    performances = evaluate_estimate(T, T_hat, Y_train2, Y_train2, Yt_train2, K)
    res_update = header.copy()
    res_update = res_update.assign(Method='Clean sample',  gamma_opt=1, **performances)
    res_list.append(res_update)

    # Loop through the methods, apply them, and evaluate the results
    for method_name, method_func in methods.items():
        print(f"Applying {method_name} method...", end=' ')
        sys.stdout.flush()

        # Initialize and apply the method for anchor points
        method = method_func()
        Ya_train2, gamma_opt, _ = method.get_anchor_points()
        if gamma_opt is None:
            gamma_opt = 50*K/n

        # Use anchor points to estimate T
        T_method = TMatrixEstimation(X_train2, Ya_train2, Yt_train2, K, estimation_method="empirical_parametricRR")
        T_hat = T_method.get_estimate()

        print("Done.")
        sys.stdout.flush()

        performances = evaluate_estimate(T, T_hat, Y_train2, Ya_train2, Yt_train2, K)
        res_update = header.copy()
        res_update = res_update.assign(Method=method_name, gamma_opt=gamma_opt, **performances)
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
