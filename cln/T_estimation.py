import numpy as np
import pandas as pd
import pdb
from scipy.spatial.distance import jensenshannon
import copy
import sys

def js_matrix_distance(P, Q):
    return np.mean([jensenshannon(P[i], Q[i]) for i in range(P.shape[0])])

def tv_matrix_distance(P, Q):
    return np.mean([0.5 * np.sum(np.abs(P[i] - Q[i])) for i in range(P.shape[0])])

def evaluate_estimate(T, T_hat, Y=None, Yt=None, K=None, anchor_points_list=None, verbose=False):

    if anchor_points_list is not None:
        # Evaluate quality of the set of anchor points
        anchors_correct = 0
        anchors_correct_tilde = 0
        anchors_total = 0

        for l in range(K):
            top_idx = anchor_points_list[l]
            Y_top = Y[top_idx]
            Yt_top = Yt[top_idx]
            anchors_correct += np.sum(Y_top == l)
            anchors_correct_tilde += np.sum(Yt_top == l)
            anchors_total   += len(top_idx)

        accuracy_gamma = anchors_correct/anchors_total
        accuracy_tilde_gamma = anchors_correct_tilde/anchors_total
    else:
        accuracy_gamma = 1
        accuracy_tilde_gamma = 1

    # Evaluate quality of estimated T
    T_norm = np.linalg.norm(T, 'fro')
    W = np.linalg.inv(T)
    Tinv_norm = np.linalg.norm(W, 'fro')

    # Method clean dataset
    #js_d = js_matrix_distance(T.T, T_hat.T)
    tv_d = tv_matrix_distance(T.T, T_hat.T)
    frobenius_d = np.linalg.norm(T - T_hat, ord='fro')/T_norm

    try:
        W_hat = np.linalg.inv(T_hat)
        frob_inv_d = np.linalg.norm(W - W_hat, ord='fro')/Tinv_norm
    except np.linalg.LinAlgError:
        frob_inv_d = None

    res = {'accuracy': accuracy_gamma,
           'accuracy_tilde': accuracy_tilde_gamma,
            'tv_d':tv_d,
            #'js_d':js_d,
            'frobenius_d':frobenius_d,
            'frob_inv_d':frob_inv_d}
         
    if verbose:
        #print('Jensen Shannon Distance:             {:2.3%}'.format(js_d))
        print('Total Variation Distance:            {:2.3%}'.format(tv_d))
        print('Fobrenius Distance:                  {:2.3%}'.format(frobenius_d))
        print('Fobrenius Distance of Inverses:      {:2.3%}'.format(frob_inv_d))

    return res

class AnchorPointsEstimation:
    def __init__(self, X, Yt, K, gamma, black_box,
                 empirical=True, verbose=False):
        
        self.K = K
        self.gamma = gamma
        self.empirical = empirical
        self.black_box = copy.deepcopy(black_box)
        self.n_cal = X.shape[0]


        # Estimate the probabilities on the calibration set
        p_hat = black_box.predict_proba(X)
        if not isinstance(p_hat, np.ndarray):
            p_hat = np.asarray(p_hat)
        _, K_out = p_hat.shape
        assert K_out == self.K, f"black_box_pt returned {K_out} classes, expected {self.K}"

        ## Estimate the contamination process by identifying anchor points
        # Identify the set of anchor points with threshold gamma%
        anchor_points_list = []

        m = max(1, int(np.ceil(self.gamma * self.n_cal)))

        T_hat = np.zeros((self.K, self.K), dtype=float)
        for l in range(self.K):
            scores_l = p_hat[:, l]
            top_idx = np.argsort(scores_l)[::-1][:m]

            # Update the list of anchor points
            anchor_points_list.append(top_idx)

            if empirical:
                # Use anchor points as true classes and then evaluate empirical frequencies
                if len(top_idx)>0:
                    counts = np.bincount(Yt[top_idx], minlength=self.K)
                    T_hat[:,l] = counts/len(top_idx)
                else:
                    T_hat[:, l] = np.ones(self.K) / self.K
            else:
                # Use the point predictor (Patrini et al.)
                col_T_hat = p_hat[top_idx, :].mean(axis=0)
                T_hat[:, l] = col_T_hat/np.sum(col_T_hat)
            
        self.T_hat = T_hat
        self.anchor_points_list = anchor_points_list

        if verbose:
            print('Estimation via anchor points completed.')
            sys.stdout.flush()

    
    def get_estimate(self):
        return self.T_hat, self.anchor_points_list