import numpy as np
from sklearn.model_selection import train_test_split
#from collections import defaultdict
from statsmodels.distributions.empirical_distribution import ECDF
import copy

from arc.classification import ProbabilityAccumulator as ProbAccum

import sys
import pdb
#import cln.optimization
from cln.optimization import estimate_c_const
from cln.optimization import eval_delta_marg
from cln.optimization import eval_delta_marg_opt

from cln.asymptotic import richardson_extrapolation, simulate_supremum

class MarginalLabelNoiseConformal:
    def __init__(self, X, Y, black_box, K, alpha, n_cal=0.5, epsilon=0.1, T=None, M=None, rho_tilde=None,
                 allow_empty=False, method="old", optimized=True, optimistic=False, 
                 asymptotic_h_start=0.0025, asymptotic_MC_samples=10000,
                 verbose=False, pre_trained=False, random_state=2023):

        self.K = K
        self.allow_empty = allow_empty
        self.method = method
        self.optimized = optimized
        self.optimistic = optimistic
        self.black_box = copy.deepcopy(black_box)
        self.asymptotic_h_start = asymptotic_h_start
        self.asymptotic_MC_samples = asymptotic_MC_samples

        if M is not None:
            self.M = M
            self.V = np.linalg.inv(M)
        
        if T is not None:
            # Use provided label noise model
            self.known_noise = True
            self.T = T
            assert T.shape[1] == K # T must be a square matrix
            self.W = np.linalg.inv(T)
            self.epsilon = epsilon
        else:
            raise ValueError("The model must be known at this stage.")

        if rho_tilde is not None:
            self.rho_tilde = rho_tilde

        # Split data into training/calibration sets
        if n_cal >=0:
            X_train, X_calib, Y_train, Y_calib = train_test_split(X, Y, test_size=n_cal, random_state=random_state)
        else:
            X_calib = X
            Y_calib = Y            
        n_cal = X_calib.shape[0]

        # Fit model (if not pre-trained)
        if pre_trained:
            if verbose:
                print('Skipping training.')
                sys.stdout.flush()
        else:
            if verbose:
                print('Fitting classifier on {:d} training samples with {:d} features...'.format(X_train.shape[0], X_train.shape[1]))
                sys.stdout.flush()
            self.black_box.fit(X_train, Y_train)
            if verbose:
                print('Training completed.')
                sys.stdout.flush()

        if verbose:
            print('Evaluating conformity scores on {:d} calibration samples...'.format(X_calib.shape[0], X_calib.shape[1]))
            sys.stdout.flush()

        # Evaluate conformity scores on calibration data
        p_hat_calib = self.black_box.predict_proba(X_calib)
        grey_box = ProbAccum(p_hat_calib)
        rng = np.random.default_rng(random_state)
        random_noise = rng.uniform(low=0.0, high=1.0, size=n_cal)
        if verbose:
            print('Evaluation of conformity scores completed.')
            sys.stdout.flush()

        if verbose:
            print('Calibrating conformity scores for {:d} classes...'.format(K))
            sys.stdout.flush()

        # Define the F-hat functions
        F_hat, scores = self.compute_F_hat_scores(grey_box, Y_calib, random_noise)
        cal_sizes = np.array([len(scores[k,0]) for k in range(self.K)])
        n_min = np.min(cal_sizes)
        rho_tilde_hat = np.divide(cal_sizes,n_cal)
        
        # Sort the conformity scores
        scores_sorted = np.concatenate([scores[(k,k)] for k in range(self.K)])
        scores_sorted = np.sort(scores_sorted)

        F_hat_marg = ECDF(scores_sorted)

        # Evaluate the Delta-hat function
        #Delta_hat = self.compute_Delta_hat_marginal(F_hat, scores_sorted)
        Delta_hat = self.compute_Delta_hat_marginal(F_hat, F_hat_marg, scores_sorted, rho_tilde_hat)
            
        # Evaluate the delta constants
        if self.method=="improved":
            delta = self.compute_delta_const_marginal_improved(n_cal)
            if not np.isscalar(delta):
                raise ValueError(f"Expected a scalar, but got {type(delta)} with value {delta}")
        elif self.method=="asymptotic":
            delta = self.compute_delta_const_marginal_asymptotic(n_cal, scores, F_hat)
        elif self.method=="old":
            delta = self.compute_delta_const_marginal(n_cal, n_min)
        else:
            raise ValueError("Unknown method for retrieving finite sample correction.")

        if self.optimistic:
            I_hat_values = np.arange(1,n_cal+1)/n_cal - (1.0-alpha) + np.maximum(-(1-alpha)/n_cal, Delta_hat - delta)
        else:
            I_hat_values = np.arange(1,n_cal+1)/n_cal - (1.0-alpha) + Delta_hat - delta
        I_hat = np.where(I_hat_values>=0)[0]

        # Calibrate the tau-hat value
        if len(I_hat)>0:
            i_hat = np.min(I_hat)
            tau_hat = scores_sorted[i_hat]
        else:
            tau_hat = 1
        tau_hat = tau_hat * np.ones((K,))

        self.tau_hat = tau_hat

        if verbose:
            print('Calibration of conformity scores completed.')
            sys.stdout.flush()

    def compute_F_hat_scores(self, grey_box, Y_calib, random_noise):
        n_cal = len(Y_calib)
        K = self.K
        scores = dict()
        F_hat = dict()
        for l in range(K):
            idx_l = np.where(Y_calib==l)[0]
            if len(idx_l)>0:
                for k in range(K):
                    # Make place-holder labels
                    Y_k = k * np.ones((n_cal,)).astype(int)
                    # Calculate conformity scores using place-holder labels
                    scores[(l,k)] = 1.0 - grey_box.calibrate_scores(Y_k, epsilon=random_noise)[idx_l]
                    F_hat[(l,k)] = ECDF(scores[(l,k)])

        return F_hat, scores

    def compute_Delta_hat_marginal(self, F_hat, F_hat_marg, scores_sorted, rho_tilde_hat):
        K = self.K
        W = self.W
        n = len(scores_sorted)
        partial_sum = np.zeros((n,))
        Delta_hat = np.zeros((n,))
        for k in range(K):
            for l in range(K):
                partial_sum += (W[k,l] * rho_tilde_hat[l]) * F_hat[(l,k)](scores_sorted)
        Delta_hat = partial_sum - F_hat_marg(scores_sorted)
        return Delta_hat

    def compute_delta_const_marginal(self, n, n_min):
        K = self.K
        c_n = estimate_c_const(n)
        tmp = np.zeros((K,))
        for k in range(K):
            for l in range(K):
                if l!=k:
                    tmp[k] += np.abs(self.V[k,l])
        rho = np.dot(self.M.T, self.rho_tilde)
        coeff = (2*np.max(tmp) + np.sum(np.abs(rho-self.rho_tilde))) / np.sqrt(n_min)
        K2 = np.power(K, 2)
        delta_k = c_n + coeff * np.minimum(K2*np.sqrt(np.pi/2), 1/np.sqrt(n_min) + np.sqrt((np.log(2*K2)+np.log(n_min))/2))
        return delta_k

    def compute_delta_const_marginal_improved(self, n):
        W = self.W
        if self.optimized is False:
            # Define the basic solution of the Randomized Response model for K
            epsilon = self.epsilon
            K = self.K
            beta_0 = 1/(1-epsilon)
            betas = np.ones((K,)) * (-epsilon/(1-epsilon))
            delta_n = eval_delta_marg(beta_0, betas, W, n)
        else:
            delta_n = eval_delta_marg_opt(W, n)
        return delta_n
    
    def compute_delta_const_marginal_asymptotic(self, n, scores, F_hat):
        W = self.W

        # Set parameters for Richardson extrapolation
        grid_type = "centered"
        h_start = self.asymptotic_h_start
        num_samples = self.asymptotic_MC_samples
        r=1

        func = lambda h: simulate_supremum(h=h, num_samples=num_samples, grid_type=grid_type,
                                           n_cal=n, scores=scores, W=W, F_hat=F_hat)
        
        exp_sup_estimate = richardson_extrapolation(h_start, func, r)

        delta_n = exp_sup_estimate/np.sqrt(n)

        return delta_n

    def predict(self, X, random_state=2023):
        n = X.shape[0]
        rng = np.random.default_rng(random_state)
        random_noise = rng.uniform(low=0.0, high=1.0, size=n)
        p_hat = self.black_box.predict_proba(X)
        grey_box = ProbAccum(p_hat)
        S_hat = [ np.array([]).astype(int) for i in range(n)]
        for k in range(self.K):
            alpha_k = 1.0 - self.tau_hat[k]
            S_k = grey_box.predict_sets(alpha_k, epsilon=random_noise, allow_empty=self.allow_empty)
            for i in range(n):
                if k in S_k[i]:
                    S_hat[i] = np.append(S_hat[i],k)
        return S_hat

