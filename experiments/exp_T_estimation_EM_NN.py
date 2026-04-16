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
from sklearn.linear_model import LogisticRegression

import sys
sys.path.append("..")
sys.path.append("../third_party")


from cln import data
from cln import contamination
from cln.T_estimation_EM import Dataset, run_em, predict
from cln.T_estimation_NN import NoisyLabelNet, train
from cln.T_estimation import evaluate_estimate


# Define default parameters
exp_num = 621
data_name = 'synthetic6'
num_var = 20
K = 4
n = 10000
n_test = 2000
n_clean = 100
model_name = 'RFC'
epsilon = 0.2
epsilon_init = 0.05
nu = 0
contamination_model = "uniform"
random_flag = True
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 12:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()

    exp_num = int(sys.argv[1])
    data_name = sys.argv[2]
    num_var = int(sys.argv[3])
    K = int(sys.argv[4])
    n = int(sys.argv[5])
    n_clean = int(sys.argv[6])
    random_flag = sys.argv[7].lower() == "true"
    epsilon = float(sys.argv[8])
    nu = float(sys.argv[9])
    contamination_model = sys.argv[10]
    seed = int(sys.argv[11])

# Define other constant parameters
batch_size = 20
n_test = 2000

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

# Initialize noise contamination process
if contamination_model == "uniform":
    T = contamination.construct_T_matrix_simple(K, epsilon)
elif contamination_model == "block":
    T = contamination.construct_T_matrix_block(K, epsilon)
elif contamination_model == "RRB":
    T = contamination.construct_T_matrix_block_RR(K, epsilon, nu)
elif contamination_model == "random":
    T = contamination.construct_T_matrix_random(K, epsilon, random_state=seed)
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

# Add important parameters to table of results
header = pd.DataFrame({'data':[data_name], 'num_var':[num_var], 'K':[K],
                       'n':[n], 'n_clean':[n_clean],
                       'epsilon':[epsilon], 'nu':[nu], 'contamination':[contamination_model],
                       'model_name':[model_name],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_p" + str(num_var)
outfile_prefix += "_K" + str(K) + "_" + model_name
outfile_prefix += "_eps" + str(epsilon) + "_epsin" + str(epsilon_init) + "_nu" + str(nu) + "_" + contamination_model
outfile_prefix += "_n" + str(n) + "_ncl" + str(n_clean) + "_seed" + str(seed)
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
    X_all, Y_all = data_distribution.sample(n + n_clean + n_test)
    print("Done.")
    sys.stdout.flush()

    # Generate the contaminated labels
    print("Generating contaminated labels...", end=' ')
    sys.stdout.flush()
    contamination_process = contamination.LinearContaminationModel(T, random_state=random_state+2)
    Yt_all = contamination_process.sample_labels(Y_all)
    print("Done.")
    sys.stdout.flush()

    # Separate data into training and calibration
    X, X_test, Y, Y_test, Yt, Yt_test = train_test_split(X_all, Y_all, Yt_all, test_size=n_test, random_state=random_state+3)

    if random_flag == True:
        rng = np.random.default_rng(random_state+4)
        clean_frac = np.round(n_clean/(n+n_clean), decimals=5)
        I = rng.binomial(1, clean_frac, size=(n+n_clean))
    else:
        n_train = 10000
        data_distribution.set_seed(random_state+5)
        X_train, Y_train = data_distribution.sample(n_train)
        Yt_train = contamination_process.sample_labels(Y_train)
        # Fit a quick RFC on noisy labels to score each observation's difficulty
        rfc_easy = RandomForestClassifier(n_estimators=100, random_state=random_state+4)
        rfc_easy.fit(X_train, Yt_train)
        conf_scores = rfc_easy.predict_proba(X).max(axis=1)

        # Assign I=1 to the top clean_frac fraction by confidence
        clean_frac = np.round(n_clean/(n+n_clean), decimals=5)
        threshold = np.quantile(conf_scores, 1 - clean_frac)
        I = (conf_scores >= threshold).astype(int)
    Y_obs = np.where(I == 1, Y, Yt)

    # Initialize an empty list to store the evaluation results
    #res = pd.DataFrame({})
    res_list = []

    ## Estimate T using the EM algorithm
    print("Estimating T using EM algorithm...", end=' ')
    sys.stdout.flush()
    X_intercept = np.hstack([np.ones((n+n_clean, 1)), X])
    data = Dataset(X=X_intercept, Y_obs=Y_obs, I=I, K=K)
    result_EM = run_em(data, eps_init=0.1, max_iter=100, tol=1e-7, verbose=False)
    eps_hat_EM = result_EM.eps
    T_hat_EM = contamination.construct_T_matrix_simple(K, eps_hat_EM)

    # predictions on test set
    X_test_intercept = np.hstack([np.ones((n_test, 1)), X_test])
    Y_test_hat_EM = predict(X_test_intercept, result_EM.beta)

    performances = evaluate_estimate(T, T_hat_EM, Y_test, Y_test_hat_EM, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='EM',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    ## Estimate T using the NN algorithm
    print("Estimating T using NN...", end=' ')
    sys.stdout.flush()
    X_torch  = torch.tensor(X, dtype=torch.float32)
    X_test_torch = torch.tensor(X_test, dtype=torch.float32)
    Y_obs_torch = torch.tensor(Y_obs, dtype=torch.long)
    I_torch = torch.tensor(I, dtype=torch.long)

    # --- Build and train model ---
    model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[32, 16], contamination_model="uniform", epsilon_init=epsilon_init)
    history = train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=1e-3, verbose=False)
    T_hat_NN = model_NN.contamination.contamination_matrix()
    T_hat_NN = T_hat_NN.detach().numpy()

    # predictions on test set
    model_NN.eval()

    with torch.no_grad():
        dummy_I     = torch.zeros(X_test_torch.shape[0], dtype=torch.long)
        dummy_noise = torch.zeros(X_test_torch.shape[0], model_NN.K)
        logits_Y, _ = model_NN(X_test_torch, dummy_I, dummy_noise)

    predicted_Y = logits_Y.argmax(dim=1)
    Y_test_hat_NN = predicted_Y.numpy()

    performances = evaluate_estimate(T, T_hat_NN, Y_test, Y_test_hat_NN, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    # Estimate T using the simplified NN algorithm
    print("Estimating T using an easier NN...", end=' ')
    sys.stdout.flush()
    model_NN_easy = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16], contamination_model="uniform", epsilon_init=epsilon_init)
    history_easy = train(model_NN_easy, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=1e-3, verbose=False)
    T_hat_NN_easy = model_NN_easy.contamination.contamination_matrix()
    T_hat_NN_easy = T_hat_NN_easy.detach().numpy()

    # predictions on test set
    model_NN_easy.eval()

    with torch.no_grad():
        dummy_I     = torch.zeros(X_test_torch.shape[0], dtype=torch.long)
        dummy_noise = torch.zeros(X_test_torch.shape[0], model_NN_easy.K)
        logits_Y_easy, _ = model_NN_easy(X_test_torch, dummy_I, dummy_noise)

    predicted_Y_easy = logits_Y_easy.argmax(dim=1)
    Y_test_hat_NN_easy = predicted_Y_easy.numpy()

    performances = evaluate_estimate(T, T_hat_NN_easy, Y_test, Y_test_hat_NN_easy, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN16',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    # Estimate T using the NN algorithm with single linear layer
    print("Estimating T using the NN with SLL...", end=' ')
    sys.stdout.flush()
    model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model="uniform", epsilon_init=epsilon_init)
    history_sll = train(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=1e-3, verbose=False)
    T_hat_NN_sll = model_NN_sll.contamination.contamination_matrix()
    T_hat_NN_sll = model_NN_sll.detach().numpy()

    # predictions on test set
    model_NN_sll.eval()

    with torch.no_grad():
        dummy_I     = torch.zeros(X_test_torch.shape[0], dtype=torch.long)
        dummy_noise = torch.zeros(X_test_torch.shape[0], model_NN_easy.K)
        logits_Y_sll, _ = model_NN_sll(X_test_torch, dummy_I, dummy_noise)

    predicted_Y_sll = logits_Y_sll.argmax(dim=1)
    Y_test_hat_NN_sll = predicted_Y_sll.numpy()

    performances = evaluate_estimate(T, T_hat_NN_sll, Y_test, Y_test_hat_NN_sll, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN SLL',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    ## Fit classifier using contaminated data and check its accuracy
    print("Performance of classifier trained on Y_obs...", end=' ')
    sys.stdout.flush()
    softmax_model = LogisticRegression(multi_class='multinomial', solver='lbfgs', max_iter=1000)
    softmax_model.fit(X, Y_obs)
    Y_test_hat_softmax = softmax_model.predict(X_test)
    performances = evaluate_estimate(T, T, Y_test, Y_test_hat_softmax, Yt_test, K, epsilon0=0)

    res_update = header.copy()
    res_update = res_update.assign(Method='softmax',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

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
