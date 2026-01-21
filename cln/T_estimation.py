import numpy as np
import pandas as pd
import pdb
import copy
import sys

def tv_matrix_distance(P, Q):
    return np.mean([0.5 * np.sum(np.abs(P[i] - Q[i])) for i in range(P.shape[0])])

def evaluate_estimate(T, T_hat, Y=None, Y_anchor_=None, Yt=None, K=None, epsilon0=0.01, verbose=False):

    # Accuracy of the set of anchor points
    accuracy = np.sum(Y_anchor_ == Y)/np.sum(Y_anchor_ != -1)
    accuracy_tilde = np.sum(Y_anchor_ == Yt)/np.sum(Y_anchor_ != -1)

    # Evaluate quality of estimated T in terms of:
    # Frobenius norm
    T_norm = np.linalg.norm(T, 'fro')
    W = np.linalg.inv(T)
    Tinv_norm = np.linalg.norm(W, 'fro')

    # Total variation distance
    tv_d = tv_matrix_distance(T.T, T_hat.T)
    frobenius_d = np.linalg.norm(T - T_hat, ord='fro')/T_norm

    # Offdiagonal mass and epsilon residual
    K = T.shape[0]
    offdiag_mass = 1/K * (np.trace(T) - np.trace(T_hat))
    epsilon_hat = (1-1/K*np.trace(T_hat))*K/(K-1)
    epsilon = (1-1/K*np.trace(T))*K/(K-1)
    epsilon_eff = epsilon + epsilon0 - epsilon*epsilon0 
    epsilon_res = epsilon_hat - epsilon_eff

    try:
        W_hat = np.linalg.inv(T_hat)
        frob_inv_d = np.linalg.norm(W - W_hat, ord='fro')/Tinv_norm
    except np.linalg.LinAlgError:
        frob_inv_d = None

    res = {'accuracy': accuracy,
           'accuracy_tilde': accuracy_tilde,
           'tv_d':tv_d,
           'frobenius_d':frobenius_d,
           'frob_inv_d':frob_inv_d,
           'offdiag_mass':offdiag_mass,
           'epsilon_res':epsilon_res}
         
    if verbose:
        print('Accuracy:             {:2.3%}'.format(accuracy))
        print('Accuracy_tilde:             {:2.3%}'.format(accuracy_tilde))
        print('Total Variation Distance:            {:2.3%}'.format(tv_d))
        print('Fobrenius Distance:                  {:2.3%}'.format(frobenius_d))
        print('Fobrenius Distance of Inverses:      {:2.3%}'.format(frob_inv_d))
        print('Off-diagonal mass:      {:2.3%}'.format(offdiag_mass))
        print('Residual with sign in epsilon estimation:      {:2.3%}'.format(epsilon_res))

    return res

class TMatrixEstimation:
    def __init__(self, X, Y_anchor, Yt, K,
                 estimation_method="empirical",
                 verbose=False):
        
        self.K = K
        self.method = estimation_method
        self.n = X.shape[0]

        # Estimate T_hat
        if estimation_method=="empirical":
            T_hat = np.zeros((self.K, self.K), dtype=float)
            for l in range(K):
                idx = (Y_anchor == l)
                n_l = np.sum(idx)
                if n_l > 0:
                    counts = np.bincount(Yt[idx], minlength=K)
                    T_hat[:, l] = counts / n_l
                else:
                    # Fallback if a class i never appears in the set of anchor points
                    T_hat[:, l] = np.ones(K) / K
            col_sums = T_hat.sum(axis=0, keepdims=True)
            T_hat /= col_sums
        # elif estimation_method=="Patrini":
        #     # Estimate the probabilities on the calibration set
        #     p_hat = black_box.predict_proba(X)
        #     if not isinstance(p_hat, np.ndarray):
        #         p_hat = np.asarray(p_hat)
        #     
        #     # Use probabilities to identify anchor points
        #     T_hat = np.zeros((self.K, self.K), dtype=float)
        #     for l in range(self.K):
        #         top_idx = self.anchor_points[l]
        #         col_T_hat = p_hat[top_idx, :].mean(axis=0)
        #         T_hat[:, l] = col_T_hat/np.sum(col_T_hat)
        elif estimation_method=="empirical_parametricRR":
            # Count how many Yt correspond to Y where Y!=-1 and use the result to estimate epsilon
            accuracy_tilde = np.sum(Y_anchor==Yt)/np.sum(Y_anchor!=-1)
            epsilon_hat = K/(K-1)*(1-accuracy_tilde)
            T_hat = (1-epsilon_hat) * np.identity(K) + epsilon_hat/K * np.ones((K,K))
            
        self.T_hat = T_hat

        if verbose:
            print('Estimation via anchor points completed.')
            sys.stdout.flush()

    def get_estimate(self):
        return self.T_hat
    
    
