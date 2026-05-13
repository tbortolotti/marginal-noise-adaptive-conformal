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
from cln.T_estimation_NN import NoisyLabelNet, train, train_alternate, train_em_style
from cln.T_estimation_EM_NN import run_em_nn
from cln.T_estimation import evaluate_estimate


# Define default parameters
exp_num = 621
data_name = 'synthetic6'
num_var = 20
K = 4
n = 10000
n_clean = 500
pi_clean = 0.5
random_flag = True
model_name = 'RFC'
epsilon = 0.2
nu = 0
contamination_model = "uniform"
contamination_exp_flag = False
seed = 1

# Parse input parameters
if True:
    print ('Number of arguments:', len(sys.argv), 'arguments.')
    print ('Argument List:', str(sys.argv))
    if len(sys.argv) != 14:
        print("Error: incorrect number of parameters.")
        quit()
    sys.stdout.flush()

    exp_num = int(sys.argv[1])
    data_name = sys.argv[2]
    num_var = int(sys.argv[3])
    K = int(sys.argv[4])
    n = int(sys.argv[5])
    n_clean = int(sys.argv[6])
    pi_clean = float(sys.argv[7])
    random_flag = sys.argv[8].lower() == "true"
    epsilon = float(sys.argv[9])
    nu = float(sys.argv[10])
    contamination_model = sys.argv[11]
    contamination_exp_flag = sys.argv[12].lower() == "true"
    seed = int(sys.argv[13])

# Define other constant parameters
batch_size = 10
epsilon_init = 0
n_test = 2000

if n_clean == 0:
    #n_clean = int(np.round(pi_clean/(1-pi_clean) * n))
    n_clean = int(np.round(pi_clean * n))
    n_noisy = n - n_clean
else:
    n_noisy = n
    n = n_noisy + n_clean

# Initialize the data distribution
if data_name == "synthetic1":
    data_distribution = data.DataModel_1(K, num_var, signal=1, random_state=seed)
elif data_name == "synthetic1_easy":
    data_distribution = data.DataModel_1_easy(K, 2, signal=1, random_state=seed)
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
                       'n':[n], 'n_noisy':[n_noisy], 'n_clean':[n_clean], 'pi_clean':[pi_clean],
                       'randflag': [random_flag],
                       'epsilon':[epsilon], 'nu':[nu], 'contamination':[contamination_model],
                       'model_name':[model_name],
                       'seed':[seed]})

# Output file
outfile_prefix = "exp"+str(exp_num) + "/" + data_name + "_p" + str(num_var)
outfile_prefix += "_K" + str(K) + "_" + model_name
outfile_prefix += "_eps" + str(epsilon) + "_epsin" + str(epsilon_init) + "_nu" + str(nu) + "_" + contamination_model
outfile_prefix += "_n" + str(n) + "_ncl" + str(n_clean) + "_picl" + str(pi_clean) + "_randf" + str(random_flag) + "_seed" + str(seed)
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
    X_all, Y_all = data_distribution.sample(n + n_test)
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
        clean_frac = np.round(n_clean/n, decimals=5)
        I = rng.binomial(1, clean_frac, size=n)
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
        clean_frac = np.round(n_clean/n, decimals=5)
        threshold = np.quantile(conf_scores, 1 - clean_frac)
        I = (conf_scores >= threshold).astype(int)
    Y_obs = np.where(I == 1, Y, Yt)

    # Initialize an empty list to store the evaluation results
    #res = pd.DataFrame({})
    res_list = []

    #____________________________________________________________________
    ## Estimate T using the EM algorithm
    print("Estimating T using EM algorithm...", end=' ')
    sys.stdout.flush()
    X_intercept = np.hstack([np.ones((n, 1)), X])
    data = Dataset(X=X_intercept, Y_obs=Y_obs, I=I, K=K)
    result_EM = run_em(data, contamination_model_="uniform", eps_init=epsilon_init, max_iter=100, tol=1e-7, verbose=False)
    T_hat_EM = result_EM.T

    # predictions on test set
    X_test_intercept = np.hstack([np.ones((n_test, 1)), X_test])
    Y_test_hat_EM = predict(X_test_intercept, result_EM.beta)

    performances = evaluate_estimate(T, T_hat_EM, Y_test, Y_test_hat_EM, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='EM',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    """
    #____________________________________________________________________
    ## Estimate T using the EM algorithm with NN
    print("Estimating T using EM algorithm with NN...", end=' ')
    sys.stdout.flush()
    data_nn = Dataset(X=X, Y_obs=Y_obs, I=I, K=K)
    result_EM_NN = run_em_nn(data_nn, contamination_model="uniform",
                            hidden_dims=[32, 16], eps_init=epsilon_init,
                            n_steps_per_iter=100, lr=1e-3, verbose=False)
    T_hat_EM_NN = result_EM_NN.T

    # predictions on test set
    Y_test_hat_EM_NN = result_EM_NN.predict(X_test)

    performances = evaluate_estimate(T, T_hat_EM_NN, Y_test, Y_test_hat_EM_NN, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN', n=n, **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    #____________________________________________________________________
    ## Estimate T using the EM algorithm with SLL
    print("Estimating T using EM algorithm with SLL...", end=' ')
    sys.stdout.flush()
    result_EM_SLL = run_em_nn(data_nn, contamination_model="uniform",
                            hidden_dims=[], eps_init=epsilon_init,
                            n_steps_per_iter=200, lr=1.0,
                            optimizer_name="lbfgs", verbose=False)
    T_hat_EM_SLL = result_EM_SLL.T

    # predictions on test set
    Y_test_hat_EM_SLL = result_EM_SLL.predict(X_test)

    performances = evaluate_estimate(T, T_hat_EM_SLL, Y_test, Y_test_hat_EM_SLL, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN SLL', n=n, **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    #____________________________________________________________________
    ## Estimate T using the EM algorithm with [16,8] NN
    print("Estimating T using EM algorithm with [16,8] NN...", end=' ')
    sys.stdout.flush()
    result_EM_SLL = run_em_nn(data_nn, contamination_model="uniform",
                            hidden_dims=[16,8], eps_init=epsilon_init,
                            n_steps_per_iter=100, lr=1e-3,
                            optimizer_name="adam", verbose=False)
    T_hat_EM_SLL = result_EM_SLL.T

    # predictions on test set
    Y_test_hat_EM_SLL = result_EM_SLL.predict(X_test)

    performances = evaluate_estimate(T, T_hat_EM_SLL, Y_test, Y_test_hat_EM_SLL, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN16', n=n, **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()
    """

    #____________________________________________________________________
    ## Estimate T using the NN algorithm
    X_torch  = torch.tensor(X, dtype=torch.float32)
    X_test_torch = torch.tensor(X_test, dtype=torch.float32)
    Y_obs_torch = torch.tensor(Y_obs, dtype=torch.long)
    I_torch = torch.tensor(I, dtype=torch.long)

    #____________________________________________________________________
    ## Estimate T using the NN with SLL and alternate training
    print("Estimating T using the NN with SLL and alt train...", end=' ')
    sys.stdout.flush()
    model_NN_sll_alt = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_sll_alt = train_alternate(model_NN_sll_alt, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=5e-2, verbose=False)
    history_sll_alt = train_alternate(model_NN_sll_alt, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

    T_hat_NN_sll_alt = model_NN_sll_alt.contamination.contamination_matrix()
    T_hat_NN_sll_alt = T_hat_NN_sll_alt.detach().numpy()

    # predictions on test set
    model_NN_sll_alt.eval()

    with torch.no_grad():
        logits_Y_sll_alt, _ = model_NN_sll_alt(X_test_torch)


    predicted_Y_sll_alt = logits_Y_sll_alt.argmax(dim=1)
    Y_test_hat_NN_sll_alt = predicted_Y_sll_alt.numpy()

    performances = evaluate_estimate(T, T_hat_NN_sll_alt, Y_test, Y_test_hat_NN_sll_alt, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN SLL alt', n=n, **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()


    #____________________________________________________________________
    ## Estimate T using the NN and alternate training
    print("Estimating T using the NN and alt train...", end=' ')
    sys.stdout.flush()
    model_NN_alt = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_alt = train_alternate(model_NN_alt, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=5e-2, verbose=False)
    history_alt = train_alternate(model_NN_alt, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)

    T_hat_NN_alt = model_NN_alt.contamination.contamination_matrix()
    T_hat_NN_alt = T_hat_NN_alt.detach().numpy()

    # predictions on test set
    model_NN_alt.eval()

    with torch.no_grad():
        logits_Y_alt, _ = model_NN_alt(X_test_torch)


    predicted_Y_alt = logits_Y_alt.argmax(dim=1)
    Y_test_hat_NN_alt = predicted_Y_alt.numpy()

    performances = evaluate_estimate(T, T_hat_NN_alt, Y_test, Y_test_hat_NN_alt, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN alt', n=n, **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    if not contamination_exp_flag:
        #____________________________________________________________________
        ## Estimate T using the NN algorithm with single linear layer
        print("Estimating T using the NN with SLL...", end=' ')
        sys.stdout.flush()
        model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
        history_sll_1 = train(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=5e-2, verbose=False)
        history_sll_2 = train(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN_sll = model_NN_sll.contamination.contamination_matrix()
        T_hat_NN_sll = T_hat_NN_sll.detach().numpy()

        # predictions on test set
        model_NN_sll.eval()

        with torch.no_grad():
            logits_Y_sll, _ = model_NN_sll(X_test_torch)

        predicted_Y_sll = logits_Y_sll.argmax(dim=1)
        Y_test_hat_NN_sll = predicted_Y_sll.numpy()

        performances = evaluate_estimate(T, T_hat_NN_sll, Y_test, Y_test_hat_NN_sll, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN SLL', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()


        #____________________________________________________________________
        ## Estimate T using the NN algorithm
        print("Estimating T using the NN...", end=' ')
        sys.stdout.flush()
        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="uniform", epsilon_init=epsilon_init)
        history_1 = train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=5e-2, verbose=False)
        history_2 = train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=1e-3, verbose=False)

        T_hat_NN = model_NN.contamination.contamination_matrix()
        T_hat_NN = T_hat_NN.detach().numpy()

        # predictions on test set
        model_NN.eval()

        with torch.no_grad():
            logits_Y, _ = model_NN(X_test_torch)

        predicted_Y = logits_Y.argmax(dim=1)
        Y_test_hat_NN = predicted_Y.numpy()

        performances = evaluate_estimate(T, T_hat_NN, Y_test, Y_test_hat_NN, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN with SLL and EM-style training
        print("Estimating T using the NN with SLL and EM-style train...", end=' ')
        sys.stdout.flush()
        model_NN_sll_ems = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
        history_sll_ems = train_em_style(model_NN_sll_ems, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128,
                                        lr_backbone=5e-2,
                                        use_closed_form_cont=True,
                                        verbose=False)
        history_sll_ems = train_em_style(model_NN_sll_ems, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128,
                                        lr_backbone=1e-3,
                                        use_closed_form_cont=True,
                                        verbose=False)

        T_hat_NN_sll_ems = model_NN_sll_ems.contamination.contamination_matrix()
        T_hat_NN_sll_ems = T_hat_NN_sll_ems.detach().numpy()

        # predictions on test set
        model_NN_sll_ems.eval()

        with torch.no_grad():
            logits_Y_sll_ems, _ = model_NN_sll_ems(X_test_torch)


        predicted_Y_sll_ems = logits_Y_sll_ems.argmax(dim=1)
        Y_test_hat_NN_sll_ems = predicted_Y_sll_ems.numpy()

        performances = evaluate_estimate(T, T_hat_NN_sll_ems, Y_test, Y_test_hat_NN_sll_ems, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN SLL ems', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()


        #____________________________________________________________________
        ## Estimate T using the NN with SLL and EM-style training
        print("Estimating T using the NN with SLL and EM-style train...", end=' ')
        sys.stdout.flush()
        model_NN_ems = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="uniform", epsilon_init=epsilon_init)
        history_ems = train_em_style(model_NN_ems, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128,
                                        lr_backbone=5e-2,
                                        use_closed_form_cont=True,
                                        verbose=False)
        history_ems = train_em_style(model_NN_ems, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128,
                                        lr_backbone=1e-3,
                                        use_closed_form_cont=True,
                                        verbose=False)

        T_hat_NN_ems = model_NN_ems.contamination.contamination_matrix()
        T_hat_NN_ems = T_hat_NN_ems.detach().numpy()

        # predictions on test set
        model_NN_ems.eval()

        with torch.no_grad():
            logits_Y_ems, _ = model_NN_ems(X_test_torch)


        predicted_Y_ems = logits_Y_ems.argmax(dim=1)
        Y_test_hat_NN_ems = predicted_Y_ems.numpy()

        performances = evaluate_estimate(T, T_hat_NN_ems, Y_test, Y_test_hat_NN_ems, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN ems', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

    """
    #____________________________________________________________________
    ## Estimate T using the NN with SLL and EM-style training
    print("Estimating T using the NN with SLL and EM-style train and lbfgs full batch...", end=' ')
    sys.stdout.flush()
    model_NN_sll_ems_lbfgs = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_sll_ems_lbfgs = train_em_style(model_NN_sll_ems_lbfgs, X_torch, Y_obs_torch, I_torch, n_epochs=200,
                                     optimizer_name="lbfgs",
                                     n_grad_steps=50, batch_size=n,
                                     lr_backbone=1.0,
                                     use_closed_form_cont=True,
                                     verbose=False)

    T_hat_NN_sll_ems_lbfgs = model_NN_sll_ems_lbfgs.contamination.contamination_matrix()
    T_hat_NN_sll_ems_lbfgs = T_hat_NN_sll_ems_lbfgs.detach().numpy()

    # predictions on test set
    model_NN_sll_ems_lbfgs.eval()

    with torch.no_grad():
        logits_Y_sll_ems_lbfgs, _ = model_NN_sll_ems_lbfgs(X_test_torch)


    predicted_Y_sll_ems_lbfgs = logits_Y_sll_ems_lbfgs.argmax(dim=1)
    Y_test_hat_NN_sll_ems_lbfgs = predicted_Y_sll_ems_lbfgs.numpy()

    performances = evaluate_estimate(T, T_hat_NN_sll_ems_lbfgs, Y_test, Y_test_hat_NN_sll_ems_lbfgs, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN SLL ems lbfgs', n=n, **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()


    #____________________________________________________________________
    ## Estimate T using the NN and EM-style training with lbfgs
    print("Estimating T using the NN and EM-style train with lbfgs...", end=' ')
    sys.stdout.flush()
    model_NN_ems_lbfgs = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_ems_lbfgs = train_em_style(model_NN_sll_ems_lbfgs, X_torch, Y_obs_torch, I_torch, n_epochs=200,
                                     optimizer_name="lbfgs",
                                     n_grad_steps=50, batch_size=n,
                                     lr_backbone=1.0,
                                     use_closed_form_cont=True,
                                     verbose=False)

    T_hat_NN_ems_lbfgs = model_NN_ems_lbfgs.contamination.contamination_matrix()
    T_hat_NN_ems_lbfgs = T_hat_NN_ems_lbfgs.detach().numpy()

    # predictions on test set
    model_NN_ems_lbfgs.eval()

    with torch.no_grad():
        logits_Y_ems_lbfgs, _ = model_NN_ems_lbfgs(X_test_torch)


    predicted_Y_ems_lbfgs = logits_Y_ems_lbfgs.argmax(dim=1)
    Y_test_hat_NN_ems_lbfgs = predicted_Y_ems_lbfgs.numpy()

    performances = evaluate_estimate(T, T_hat_NN_ems_lbfgs, Y_test, Y_test_hat_NN_ems_lbfgs, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NN ems lbfgs', n=n, **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()
    """

    if contamination_exp_flag:
        #____________________________________________________________________
        ## Estimate T using the EM algorithm with general contamination model
        print("Estimating T using EM algorithm with general contamination...", end=' ')
        sys.stdout.flush()
        X_intercept = np.hstack([np.ones((n, 1)), X])
        data = Dataset(X=X_intercept, Y_obs=Y_obs, I=I, K=K)
        result_EM = run_em(data, contamination_model_="general", eps_init=epsilon_init, max_iter=100, tol=1e-7, verbose=False)
        T_hat_EM = result_EM.T

        # predictions on test set
        X_test_intercept = np.hstack([np.ones((n_test, 1)), X_test])
        Y_test_hat_EM = predict(X_test_intercept, result_EM.beta)

        performances = evaluate_estimate(T, T_hat_EM, Y_test, Y_test_hat_EM, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='EM gen',  **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the EM algorithm with NN and general contamination model
        print("Estimating T using EM algorithm with NN and general contamination...", end=' ')
        sys.stdout.flush()
        data_nn = Dataset(X=X, Y_obs=Y_obs, I=I, K=K)
        result_EM_NN = run_em_nn(data_nn, contamination_model="general",
                            hidden_dims=[16,8], eps_init=epsilon_init,
                            n_steps_per_iter=100,
                            lr=1e-3,
                            optimizer_name="adam",
                            verbose=False)
        T_hat_EM_NN = result_EM_NN.T

        # predictions on test set
        Y_test_hat_EM_NN = result_EM_NN.predict(X_test)

        performances = evaluate_estimate(T, T_hat_EM_NN, Y_test, Y_test_hat_EM_NN, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='EM NN gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN algorithm with general contamination
        print("Estimating T using NN with general contamination...", end=' ')
        sys.stdout.flush()
        X_torch  = torch.tensor(X, dtype=torch.float32)
        X_test_torch = torch.tensor(X_test, dtype=torch.float32)
        Y_obs_torch = torch.tensor(Y_obs, dtype=torch.long)
        I_torch = torch.tensor(I, dtype=torch.long)

        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16, 8], contamination_model_="general", epsilon_init=epsilon_init)
        train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=5e-2, verbose=False)
        train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN = model_NN.contamination.contamination_matrix()
        T_hat_NN = T_hat_NN.detach().numpy()

        # predictions on test set
        model_NN.eval()

        with torch.no_grad():
            logits_Y, _ = model_NN(X_test_torch)

        predicted_Y = logits_Y.argmax(dim=1)
        Y_test_hat_NN = predicted_Y.numpy()

        performances = evaluate_estimate(T, T_hat_NN, Y_test, Y_test_hat_NN, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN algorithm with alt training and general contamination
        print("Estimating T using NN with alt training and general contamination...", end=' ')
        sys.stdout.flush()

        model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16, 8], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=5e-2, verbose=False)
        train_alternate(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN = model_NN.contamination.contamination_matrix()
        T_hat_NN = T_hat_NN.detach().numpy()

        # predictions on test set
        model_NN.eval()

        with torch.no_grad():
            logits_Y, _ = model_NN(X_test_torch)

        predicted_Y = logits_Y.argmax(dim=1)
        Y_test_hat_NN = predicted_Y.numpy()

        performances = evaluate_estimate(T, T_hat_NN, Y_test, Y_test_hat_NN, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN alt gen',  **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN algorithm with single linear layer with general contamination
        print("Estimating T using the NN with SLL and general contamination...", end=' ')
        sys.stdout.flush()
        model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="general", epsilon_init=epsilon_init)
        train_alternate(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=5e-2, verbose=False)
        train_alternate(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=100, n_grad_steps=50, batch_size=128, lr=1e-3, verbose=False)
        T_hat_NN_sll = model_NN_sll.contamination.contamination_matrix()
        T_hat_NN_sll = T_hat_NN_sll.detach().numpy()

        # predictions on test set
        model_NN_sll.eval()

        with torch.no_grad():
            logits_Y_sll, _ = model_NN_sll(X_test_torch)

        predicted_Y_sll = logits_Y_sll.argmax(dim=1)
        Y_test_hat_NN_sll = predicted_Y_sll.numpy()

        performances = evaluate_estimate(T, T_hat_NN_sll, Y_test, Y_test_hat_NN_sll, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN SLL alt gen',  **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()

        #____________________________________________________________________
        ## Estimate T using the NN and EM-style training with general contamination
        print("Estimating T using the NN and EM-style train with general contamination...", end=' ')
        sys.stdout.flush()
        model_NN_ems = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[16,8], contamination_model_="general", epsilon_init=epsilon_init)
        train_em_style(model_NN_ems, X_torch, Y_obs_torch, I_torch, n_epochs=100, batch_size=128, lr_backbone=1e-3, verbose=False)

        T_hat_NN_ems = model_NN_ems.contamination.contamination_matrix()
        T_hat_NN_ems = T_hat_NN_ems.detach().numpy()

        # predictions on test set
        model_NN_ems.eval()

        with torch.no_grad():
            logits_Y_ems, _ = model_NN_ems(X_test_torch)

        predicted_Y_ems = logits_Y_ems.argmax(dim=1)
        Y_test_hat_NN_ems = predicted_Y_ems.numpy()

        performances = evaluate_estimate(T, T_hat_NN_ems, Y_test, Y_test_hat_NN_ems, Yt_test, K, epsilon0=0)
        res_update = header.copy()
        res_update = res_update.assign(Method='NN ems gen', n=n, **performances)
        res_list.append(res_update)
        print("Done.")
        sys.stdout.flush()


    """
    #____________________________________________________________________
    ## Estimate T using the NN algorithm with weighted loss
    print("Estimating T using NN with weighted loss...", end=' ')
    sys.stdout.flush()
    model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[32, 16], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_w_1 = train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=5e-2, loss_type="weighted", verbose=False)
    history_w_2 = train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=1e-3, loss_type="weighted", verbose=False)
    T_hat_NNw = model_NN.contamination.contamination_matrix()
    T_hat_NNw = T_hat_NNw.detach().numpy()

    # predictions on test set
    model_NN.eval()

    with torch.no_grad():
        dummy_I     = torch.zeros(X_test_torch.shape[0], dtype=torch.long)
        dummy_noise = torch.zeros(X_test_torch.shape[0], model_NN.K)
        logits_w_Y, _ = model_NN(X_test_torch, dummy_I, dummy_noise)

    predicted_w_Y = logits_w_Y.argmax(dim=1)
    Y_test_hat_NNw = predicted_w_Y.numpy()

    performances = evaluate_estimate(T, T_hat_NNw, Y_test, Y_test_hat_NNw, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NNw',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    #____________________________________________________________________
    ## Estimate T using the NN algorithm with single linear layer with weighted loss
    print("Estimating T using the NN with SLL and weighted loss...", end=' ')
    sys.stdout.flush()
    model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_sllw_1 = train(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=5e-2, loss_type="weighted", verbose=False)
    history_sllw_2 = train(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=1e-3, loss_type="weighted", verbose=False)
    T_hat_NN_sllw = model_NN_sll.contamination.contamination_matrix()
    T_hat_NN_sllw = T_hat_NN_sllw.detach().numpy()

    # predictions on test set
    model_NN_sll.eval()

    with torch.no_grad():
        dummy_I     = torch.zeros(X_test_torch.shape[0], dtype=torch.long)
        dummy_noise = torch.zeros(X_test_torch.shape[0], model_NN_sll.K)
        logits_Y_sllw, _ = model_NN_sll(X_test_torch, dummy_I, dummy_noise)

    predicted_Y_sllw = logits_Y_sllw.argmax(dim=1)
    Y_test_hat_NN_sllw = predicted_Y_sllw.numpy()

    performances = evaluate_estimate(T, T_hat_NN_sllw, Y_test, Y_test_hat_NN_sllw, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NNw SLL',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    #____________________________________________________________________
    ## Estimate T using the NN algorithm with upweighted loss
    print("Estimating T using NN with upweighted loss...", end=' ')
    sys.stdout.flush()
    model_NN = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[32, 16], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_uw_1 = train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=1e-2, loss_type="upweighted", verbose=False)
    history_uw_2 = train(model_NN, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=1e-3, loss_type="upweighted", verbose=False)
    T_hat_NNuw = model_NN.contamination.contamination_matrix()
    T_hat_NNuw = T_hat_NNuw.detach().numpy()

    # predictions on test set
    model_NN.eval()

    with torch.no_grad():
        dummy_I     = torch.zeros(X_test_torch.shape[0], dtype=torch.long)
        dummy_noise = torch.zeros(X_test_torch.shape[0], model_NN.K)
        logits_uw_Y, _ = model_NN(X_test_torch, dummy_I, dummy_noise)

    predicted_uw_Y = logits_uw_Y.argmax(dim=1)
    Y_test_hat_NNuw = predicted_uw_Y.numpy()

    performances = evaluate_estimate(T, T_hat_NNuw, Y_test, Y_test_hat_NNuw, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NNuw',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()

    #____________________________________________________________________
    ## Estimate T using the NN algorithm with single linear layer with upweighted loss
    print("Estimating T using the NN with SLL and upweighted loss...", end=' ')
    sys.stdout.flush()
    model_NN_sll = NoisyLabelNet(input_dim=num_var, K=K, hidden_dims=[], contamination_model_="uniform", epsilon_init=epsilon_init)
    history_slluw_1 = train(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=1e-2, loss_type="upweighted", verbose=False)
    history_slluw_2 = train(model_NN_sll, X_torch, Y_obs_torch, I_torch, n_epochs=500, batch_size=128, lr=1e-3, loss_type="upweighted", verbose=False)
    T_hat_NN_slluw = model_NN_sll.contamination.contamination_matrix()
    T_hat_NN_slluw = T_hat_NN_slluw.detach().numpy()

    # predictions on test set
    model_NN_sll.eval()

    with torch.no_grad():
        dummy_I     = torch.zeros(X_test_torch.shape[0], dtype=torch.long)
        dummy_noise = torch.zeros(X_test_torch.shape[0], model_NN_sll.K)
        logits_Y_slluw, _ = model_NN_sll(X_test_torch, dummy_I, dummy_noise)

    predicted_Y_slluw = logits_Y_slluw.argmax(dim=1)
    Y_test_hat_NN_slluw = predicted_Y_slluw.numpy()

    performances = evaluate_estimate(T, T_hat_NN_slluw, Y_test, Y_test_hat_NN_slluw, Yt_test, K, epsilon0=0)
    res_update = header.copy()
    res_update = res_update.assign(Method='NNuw SLL',  **performances)
    res_list.append(res_update)
    print("Done.")
    sys.stdout.flush()
    """

    #____________________________________________________________________
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
