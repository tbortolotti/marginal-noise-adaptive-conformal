import numpy as np

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
from cln.AP_identification import AnchorPointsExistence
from third_party import arc

# Define default parameters
exp_num = 611

# Data simulation parameters
data_name = 'syntheticAP'
scenario = "scenario1"

# Contamination parameters
contamination_model = "uniform"
epsilon = 0.1

# Set sizes
n_train1 = 10000
n_train2 = 5000
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
    data_name = sys.argv[2]
    scenario = sys.argv[3]
    epsilon = float(sys.argv[4])
    contamination_model = sys.argv[5]
    n_train1 = int(sys.argv[6])
    n_train2 = int(sys.argv[7])
    seed = int(sys.argv[8])

# Define other parameters
K = 4
num_var = 2
sigma_easy = 1
sigma_hard = 1
R_easy = 1.5
R_hard = 4.0
flipy = 0
A_easy = sigma_easy * np.eye(2)
A_hard = sigma_hard * np.eye(2)
nu = 0

if scenario == "scenario1":
    # Scenario 1: Anchor points exist
    pi_easy=1
    delta_shift=0
    center_scale=0.75
elif scenario == "scenario2":
    # Scenario 2: Anchor points do not exist
    pi_easy=0.5
    delta_shift=0
    center_scale=0.75
elif scenario == "scenario3":
    # Scenario 3: Absence of anchor points is more marked than in Scenario3
    pi_easy=0.5
    delta_shift=0.25
    center_scale=0.75

batch_size = 10

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
    print("Unknown contamination model!")
    sys.stdout.flush()
    exit(-1)

# Initialize black-box model
black_box_SVC = arc.black_boxes.SVC(clip_proba_factor = 1e-5)

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'scenario':[scenario],
                       'n_train1':[n_train1], 'n_train2':[n_train2],
                       'epsilon':[epsilon], 'contamination':[contamination_model],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_" + scenario
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

    methods = {

        "Split 5": lambda: AnchorPointsExistence(X, Yt, K, n2=n_train2, black_box=black_box_SVC,
                                             method="split", frac_threshold = 0.05, consistency_threshold=0.6),

        #"Split 10": lambda: AnchorPointsExistence(X, Yt, K, n2=n_train2, black_box=black_box_SVC,
        #                                     method="split", frac_threshold = 0.1, consistency_threshold=0.6),

        #"Split 20": lambda: AnchorPointsExistence(X, Yt, K, n2=n_train2, black_box=black_box_SVC,
        #                                     method="split", frac_threshold = 0.2, consistency_threshold=0.6),

        "Boot 5": lambda: AnchorPointsExistence(X, Yt, K, n2=n_train2, black_box=black_box_SVC,
                                             method="bootstrap", freq_threshold = 0.3, stable_frac_threshold=0.05)

        #"Boot 10": lambda: AnchorPointsExistence(X, Yt, K, n2=n_train2, black_box=black_box_SVC,
        #                                     method="bootstrap", freq_threshold = 0.3, stable_frac_threshold=0.1),

        #"Boot 20": lambda: AnchorPointsExistence(X, Yt, K, n2=n_train2, black_box=black_box_SVC,
        #                                     method="bootstrap", freq_threshold = 0.3, stable_frac_threshold=0.2),
    }

    # Initialize an empty list to store the evaluation results
    res = pd.DataFrame({})
    res_list = []

    # Loop through the methods, apply them, and evaluate the results
    for method_name, method_func in methods.items():
        print(f"Applying {method_name} method...", end=' ')
        sys.stdout.flush()

        # Initialize and apply the method for assessing the existence of anchor points
        method = method_func()
        existence = method.get_anchors_exist()

        print("Done.")
        sys.stdout.flush()

        res_update = header.copy()
        res_update = res_update.assign(
            Method       = method_name,
            n_train1     = n_train1,
            n_train2     = n_train2, 
            existence    = existence,
            truth        = (scenario == "scenario1"),
            correct      = existence == (scenario == "scenario1"),
            FP           = existence and (scenario != "scenario1"),
            FN           = (not existence) and (scenario == "scenario1"),
        )
        res_update = res_update.assign(Method=method_name, existence=existence)
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
