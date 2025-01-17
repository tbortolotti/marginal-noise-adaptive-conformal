import numpy as np
from sklearn.model_selection import train_test_split

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm
import pdb

import sys
sys.path.append("..")
sys.path.append("../third_party")


from cln import data
from cln import contamination
from cln.utils import evaluate_predictions, estimate_rho

from cln.classification import MarginalLabelNoiseConformal
from cln.classification_label_conditional import LabelNoiseConformal

from third_party import arc


# Define default parameters
exp_num = 1
data_name = 'synthetic4'
num_var = 20
K = 4
signal = 1
model_name = 'RFC'
epsilon = 0.1
nu = 0
contamination_model = "uniform"
n_train = 1000
n_cal = 5000
estimate = "None"
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 15:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()

    exp_num = int(sys.argv[1])
    data_name = sys.argv[2]
    num_var = int(sys.argv[3])
    K = int(sys.argv[4])
    signal = float(sys.argv[5])
    model_name = sys.argv[6]
    epsilon = float(sys.argv[7])
    nu = float(sys.argv[8])
    contamination_model = sys.argv[9]
    n_train = int(sys.argv[10])
    n_cal = int(sys.argv[11])
    estimate = sys.argv[12]
    imb = int(sys.argv[13])
    seed = int(sys.argv[14])


# Define other constant parameters
n_test = 2000
batch_size = 20
allow_empty = True
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000

# Initialize the data distribution
if data_name == "synthetic1":
    data_distribution = data.DataModel_1(K, num_var, signal=signal, random_state=seed)
elif data_name == "synthetic2":
    data_distribution = data.DataModel_2(K, num_var, signal=signal, random_state=seed)
elif data_name == "synthetic3":
    data_distribution = data.DataModel_3(K, num_var, signal=signal, random_state=seed)
elif data_name == "synthetic4":
    data_distribution = data.DataModel_4(K, num_var, imb=imb, random_state=seed)
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


# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'num_var':[num_var], 'K':[K],
                       'signal':[signal], 'n_train':[n_train], 'n_cal':[n_cal],
                       'epsilon':[epsilon], 'nu':[nu], 'contamination':[contamination_model],
                       'model_name':[model_name], 'estimate':[estimate], 'imb':[imb], 'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_p" + str(num_var)
outfile_prefix += "_K" + str(K) + "_signal" + str(signal) + "_" + model_name
outfile_prefix += "_eps" + str(epsilon) + "_nu" + str(nu) + "_" + contamination_model
outfile_prefix += "_nt" + str(n_train) + "_nc" + str(n_cal) + "_est" + estimate + "_imb" + str(imb) +  "_seed" + str(seed)
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
        rho_hat = rho
        rho_tilde_hat = rho_tilde
        M_hat = M
    elif estimate=="rho":
        rho_tilde_hat = estimate_rho(Yt, K)
        rho_hat = np.dot(M.T, rho_tilde_hat)
        M_hat = M
    else:
        print("Unknown estimation option!")
        sys.stdout.flush()
        exit(-1)


    # Apply standard method to corrupted labels (for training)
    print("Training the predictive model...", end=' ')
    sys.stdout.flush()
    method_train = arc.methods.SplitConformal(X, Yt, black_box, K, 0.1, n_cal=n_cal, random_state=random_state)
    print("Done.")
    sys.stdout.flush()

    # Extract the pre-trained model
    black_box_pt = method_train.black_box

    res = pd.DataFrame({})
    for alpha in [0.1]:
        #for guarantee in ['lab-cond', 'marginal']:
        for guarantee in ['marginal']:
            print("\nSeeking {:s} coverage at level {:.2f}.".format(guarantee, 1-alpha))

            #if guarantee=='lab-cond':
            #    label_conditional = True
            #else:
            label_conditional = False

            # Define a dictionary of methods with their names and corresponding initialization parameters
            methods = {
                "Standard": lambda: arc.methods.SplitConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                               label_conditional=label_conditional, allow_empty=allow_empty,
                                                               pre_trained=True, random_state=random_state),

                "Adaptive optimized": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                          epsilon=epsilon, T=T, rho_tilde=rho_tilde_hat,
                                                                          allow_empty=allow_empty, method="improved",
                                                                          optimized=True, optimistic=False, verbose=False,
                                                                          pre_trained=True, random_state=random_state),

                "Adaptive optimized+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                          epsilon=epsilon, T=T, rho_tilde=rho_tilde_hat,
                                                                          allow_empty=allow_empty, method="improved",
                                                                          optimized=True, optimistic=True, verbose=False,
                                                                          pre_trained=True, random_state=random_state),

                "Adaptive simplified": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                           epsilon=epsilon, T=T, rho_tilde=rho_tilde_hat,
                                                                           allow_empty=allow_empty, method="improved",
                                                                           optimized=False, optimistic=False, verbose=False,
                                                                           pre_trained=True, random_state=random_state),

                "Adaptive simplified+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                           epsilon=epsilon, T=T, rho_tilde=rho_tilde_hat,
                                                                           allow_empty=allow_empty, method="improved",
                                                                           optimized=False, optimistic=True, verbose=False,
                                                                           pre_trained=True, random_state=random_state),

                "Asymptotic": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                  epsilon=epsilon, asymptotic_h_start=asymptotic_h_start,
                                                                  asymptotic_MC_samples=asymptotic_MC_samples, T=T,
                                                                  rho_tilde=rho_tilde_hat, allow_empty=allow_empty,
                                                                  method="asymptotic", optimistic=False, verbose=False,
                                                                  pre_trained=True, random_state=random_state),

                "Asymptotic+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                   epsilon=epsilon, asymptotic_h_start=asymptotic_h_start,
                                                                   asymptotic_MC_samples=asymptotic_MC_samples, T=T,
                                                                   rho_tilde=rho_tilde_hat, allow_empty=allow_empty,
                                                                   method="asymptotic", optimistic=True, verbose=False,
                                                                   pre_trained=True, random_state=random_state),

                "Standard label conditional": lambda: arc.methods.SplitConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                               label_conditional=True, allow_empty=allow_empty,
                                                               pre_trained=True, random_state=random_state),

                "Label conditional": lambda: LabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                 rho_tilde=rho_tilde_hat, M=M_hat,
                                                                 calibration_conditional=False, gamma=None,
                                                                 optimistic=False, allow_empty=allow_empty, verbose=False, pre_trained=True, random_state=random_state),
                
                "Label conditional+": lambda: LabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                  rho_tilde=rho_tilde_hat, M=M_hat,
                                                                  calibration_conditional=False, gamma=None,
                                                                  optimistic=True, allow_empty=allow_empty, verbose=False, pre_trained=True, random_state=random_state)

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
