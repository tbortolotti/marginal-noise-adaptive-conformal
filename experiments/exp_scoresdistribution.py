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
data_name = 'synthetic1'
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
    if len(sys.argv) != 13:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()

    exp_num = int(sys.argv[1])
    data_name = sys.argv[2]
    num_var = int(sys.argv[3])
    K = int(sys.argv[4])
    model_name = sys.argv[5]
    epsilon = float(sys.argv[6])
    nu = float(sys.argv[7])
    contamination_model = sys.argv[8]
    n_train = int(sys.argv[9])
    n_cal = int(sys.argv[10])
    estimate = sys.argv[11]
    seed = int(sys.argv[12])


# Define other constant parameters
n_test = 2000
batch_size = 20
allow_empty = True
asymptotic_h_start = 1/400
asymptotic_MC_samples = 10000

m = 2
lambda_ = 0.8
var_ = 0.02
geometry = "circle"

def sample_by_label_distribution(X_all, Y_all, rho, total_samples, random_state = 2025):
    """
    X_all:          numpy array of shape (n, p)
    Y_all:          list or 1D array of length n (categorical labels)
    rho:            dict or mapping {label: fraction}, e.g., {0: 0.2, 1: 0.5, 2: 0.3}
    total_samples:  integer, the size of the final sampled dataset

    Returns:
    X_sampled, Y_sampled: the new subset of data
    """
    # Convert Y_all to a NumPy array if it's not already
    Y_all = np.array(Y_all)

    # Prepare lists to hold sampled X and Y
    X_subsets = []
    Y_subsets = []
    rng = np.random.default_rng(seed=random_state)

    # For each label, sample according to its fraction in rho
    for label, fraction in rho.items():
        # Indices in Y_all for this label
        label_indices = np.where(Y_all == int(label))[0]

        # Number of samples needed for this label
        num_samples_for_label = int(round(fraction * total_samples))
        # pdb.set_trace()

        # Edge case: if label_indices is empty or fraction is 0
        if len(label_indices) == 0 or num_samples_for_label == 0:
            continue

        # Randomly choose the required number of indices
        chosen_indices = rng.choice(label_indices,
                                          size=num_samples_for_label,
                                          replace=False)


        # Gather the corresponding rows from X_all and labels from Y_all
        X_subsets.append(X_all[chosen_indices, :])
        Y_subsets.append(Y_all[chosen_indices])

    # Concatenate everything
    X_sampled = np.concatenate(X_subsets, axis=0)
    Y_sampled = np.concatenate(Y_subsets, axis=0)

    return X_sampled, Y_sampled

def generate_synthetic_data(p, K, m, lambda_,
                            n_each=100,      # points per cluster
                            var=0.05,    # cluster variance
                            spacing=1,# spacing between original cluster centers
                            geometry = "line",
                            random_state=2025,
                            randomize = True,
                            verbose = True,
                            return_num = False):
    """
    Generate synthetic data in a p-dimensional space with K clusters.

    Rules:
      1. Each cluster has variance = var (isotropic in p-dimensions).
      2. When lambda_ = 0, all K cluster centers are equally spaced.
      3. As lambda_ increases (0 <= lambda_ <= 1), randomly choose m clusters
         as 'main clusters'. The remaining (K - m) clusters are pulled toward
         their nearest main cluster by a factor of lambda_.

    Parameters
    ----------
    p : int
        Number of features (dimensions).
    K : int
        Number of clusters (classes).
    m : int
        Number of main clusters among the K clusters.
    lambda_ : float
        The parameter controlling how much the non-main clusters move closer
        to their nearest main cluster. (0 = no movement, 1 = fully merged)
    n : int, optional (default=100)
        Number of points per cluster.
    var : float, optional (default=1.0)
        Variance for the isotropic Gaussian that generates each cluster.
    spacing : float, optional (default=5.0)
        The initial distance between consecutive cluster centers.
    geometry : str, optional (default = "circle)
        The initial geometry of the cluster centers.
    random_state : int or None, optional (default=42)
        Random seed for reproducibility. If None, no fixed seed is used.

    Returns
    -------
    X : np.ndarray of shape (K*n, p)
        The generated feature matrix.
    Y : np.ndarray of shape (K*n,)
        The integer cluster labels for each point.
    """
    assert m <= K, "The number of main clusters (m) must be smaller than or equal to the number of total clusters (K)"
    assert 0 <= lambda_ <=1, "The similarity parameter lambda must be between 0 and 1"

    rng = np.random.default_rng(random_state)
    centers = []

    if geometry == "line":# initial cluster centers: [ i * spacing, 0, 0,..0] in p-dimensional space.
      for i in range(K):
          center = np.zeros(p, dtype=float)
          center[0] = i * spacing
          centers.append(center)


    elif geometry == "circle":
      offsets = rng.uniform(-spacing*var, spacing*var, size=K)
      angles = np.linspace(0, 2*np.pi, K, endpoint=False)
      for i in range(K):
          center = np.zeros(p, dtype=float)
          center[0] = spacing * np.cos(angles[i]) + offsets[i]
          center[1] = spacing * np.sin(angles[i]) + offsets[i]
          centers.append(center)

    # elif geometry == "hypercube":

    centers = np.array(centers)

    # Randomly choose m "main clusters"
    # Then for each of the remaining (K-m) clusters, move them closer to the nearest main cluster by a factor of frac_i
    # final_centers[i] = (1 - frac_i) * final_centers[i] + frac_i * main_centers[nearest_main_idx]

    def find_nearest_main(center_i):
        dists_sq = np.sum((main_centers - center_i)**2, axis=1)
        return np.argmin(dists_sq)

    final_centers = centers.copy()

    offsets = rng.uniform(0, 0.1, size=K)
    main_idx = [2, 7]# rng.choice(K, size=m, replace=False)
    main_centers = centers[main_idx]

    if lambda_ > 0:

        for i in range(K):
            if i not in main_idx:
                nearest_main_idx = find_nearest_main(centers[i])

                frac_i = min(lambda_ * (1 + offsets[i]), 0.99) if randomize else lambda_
                final_centers[i] = (1 - frac_i) * final_centers[i] + frac_i * main_centers[nearest_main_idx]

    # Generate data points for each cluster from N(final_center_i, var * I)
    X = np.zeros((K*n_each, p))
    Y = np.zeros(K*n_each, dtype=int)

    start_idx = 0
    for i in range(K):
        cov = var * np.eye(p)
        points = rng.multivariate_normal(mean = final_centers[i], cov = cov, size = n_each)
        X[start_idx : start_idx + n_each] = points
        Y[start_idx : start_idx + n_each] = i
        start_idx += n_each

    if verbose:
      if lambda_ > 0:
        print(f"Geometry: {geometry} \n {m} number of main clusters with centers: \n {main_centers}")
      else:
        print(f"Geometry: {geometry}")

    return (X, Y, final_centers, main_idx) if return_num else (X, Y)

# Impose the label proportions from the population model
rho_list = {
    "0": 0.02,
    "1": 0.15,
    "2": 0.13,
    "3": 0.09,
    "4": 0.06,
    "5": 0.05,
    "6": 0.04,
    "7": 0.03,
    "8": 0.4,
    "9": 0.03
}

rho = [0.02, 0.15, 0.13, 0.09, 0.06, 0.05, 0.04, 0.03, 0.4, 0.03]

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
                       'model_name':[model_name], 'estimate':[estimate], 'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_p" + str(num_var)
outfile_prefix += "_K" + str(K) + "_signal" + str(signal) + "_" + model_name
outfile_prefix += "_eps" + str(epsilon) + "_nu" + str(nu) + "_" + contamination_model
outfile_prefix += "_nt" + str(n_train) + "_nc" + str(n_cal) + "_est" + estimate + "_seed" + str(seed)
print("Output file: {:s}.".format("results/"+outfile_prefix), end="\n")
sys.stdout.flush()

# Describe the experiment
def run_experiment(random_state):
    print("\nRunning experiment in batch {:d}...".format(random_state))
    sys.stdout.flush()

    # Generate a large data set
    print("\nGenerating data...", end=' ')
    sys.stdout.flush()
    X_all, Y_all = generate_synthetic_data(p=num_var, K=K, m=m, lambda_= lambda_,
                            n_each=n_train + n_cal + n_test,
                            var=var_,
                            spacing=1,
                            geometry = geometry,
                            random_state=random_state,
                            randomize = True)
    X_all, Y_all = sample_by_label_distribution(X_all, Y_all, rho_list, total_samples = n_train+n_cal+n_test,random_state = random_state+1)
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

                "Adaptive": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                epsilon=epsilon, T=T, M=M, rho_tilde=rho_tilde_hat,
                                                                allow_empty=allow_empty, method="old", optimistic=False,
                                                                verbose=False, pre_trained=True, random_state=random_state),

                "Adaptive+": lambda: MarginalLabelNoiseConformal(X, Yt, black_box_pt, K, alpha, n_cal=n_cal,
                                                                epsilon=epsilon, T=T, M=M, rho_tilde=rho_tilde_hat,
                                                                allow_empty=allow_empty, method="old", optimistic=True,
                                                                verbose=False, pre_trained=True, random_state=random_state),

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
