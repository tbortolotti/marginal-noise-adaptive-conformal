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
    K = T.shape[0]
    offdiag_mass = 1/K * (np.trace(T) - np.trace(T_hat))

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
            'frob_inv_d':frob_inv_d,
            'offdiag_mass':offdiag_mass}
         
    if verbose:
        #print('Jensen Shannon Distance:             {:2.3%}'.format(js_d))
        print('Total Variation Distance:            {:2.3%}'.format(tv_d))
        print('Fobrenius Distance:                  {:2.3%}'.format(frobenius_d))
        print('Fobrenius Distance of Inverses:      {:2.3%}'.format(frob_inv_d))

    return res

def evaluate_accuracy_tilde(Yt, K, anchor_points_):
    # Evaluate quality of the set of anchor points
    anchors_correct_tilde = 0
    anchors_total = 0

    for l in range(K):
        top_idx = anchor_points_[l]
        Yt_top = Yt[top_idx]
        anchors_correct_tilde += np.sum(Yt_top == l)
        anchors_total   += len(top_idx)
    
    accuracy_tilde_gamma = anchors_correct_tilde/anchors_total

    return accuracy_tilde_gamma

def gamma_by_percent_drop(gamma_vec, accuracy_vec, drop=0.03):
    acc_max = np.max(accuracy_vec)
    ok = np.where(accuracy_vec >= (1.0-drop)*acc_max)[0]
    if ok.size == 0:
        idx = int(np.argmax(accuracy_vec))
    else:
        idx = int(ok[-1])
    return gamma_vec[idx]

"""
def elbow_gamma_distance_to_line(gamma_vec, accuracy_tilde_vec, smooth_flag=False):
    # Distance-to-line elbow finder.
    if gamma_vec.size != accuracy_tilde_vec.size or gamma_vec.size < 3:
        raise ValueError("gamma_vec and accuracy_tilde_vec must have same length >= 3.")
    
    y = accuracy_tilde_vec

    if smooth_flag:
        window = 3
        pad = window // 2
        ypad = np.pad(y, (pad, pad), mode="edge")
        kernel = np.ones(window) / window
        y = np.convolve(ypad, kernel, mode="valid")

    # Endpoints
    x1, y1 = gamma_vec[0], y[0]
    x2, y2 = gamma_vec[-1], y[-1]

    dx, dy = (x2 - x1), (y2 - y1)
    denom = np.hypot(dx, dy)
    if denom == 0:
        raise ValueError("Endpoints are identical; cannot form a line for distance-to-line method.")

    # Perpendicular distances
    distances = np.abs((y1 - y2) * gamma_vec + (x2 - x1) * y + (x1 * y2 - x2 * y1)) / denom

    # Exclude endpoints
    distances[0] = -np.inf
    distances[-1] = -np.inf

    idx_elbow = int(np.argmax(distances))
    gamma_elbow = gamma_vec[idx_elbow]

    return gamma_elbow
"""

def elbow_gamma_distance_to_line(gamma_vec, accuracy_tilde_vec, smooth_flag=False):
    # Distance-to-line elbow finder in log(gamma) space.

    # Convert to numeric arrays
    gamma_vec = np.asarray(gamma_vec, dtype=float).ravel()
    y = np.asarray(accuracy_tilde_vec, dtype=float).ravel()

    if smooth_flag:
        window = 3
        pad = window // 2
        ypad = np.pad(y, (pad, pad), mode="edge")
        kernel = np.ones(window) / window
        y = np.convolve(ypad, kernel, mode="valid")

    if gamma_vec.size != y.size or gamma_vec.size < 3:
        raise ValueError("gamma_vec and accuracy_tilde_vec must have same length >= 3.")

    if np.any(gamma_vec <= 0):
        raise ValueError("All gamma values must be strictly positive to take logs.")

    # Log-transform gamma
    x = np.log10(gamma_vec)

    # Endpoints of the chord
    x1, y1 = x[0], y[0]
    x2, y2 = x[-1], y[-1]

    dx, dy = (x2 - x1), (y2 - y1)
    denom = np.hypot(dx, dy)
    if denom == 0:
        raise ValueError("Endpoints are identical; cannot form a line for distance-to-line method.")

    # Line coefficients: A x + B y + C = 0
    A = y1 - y2
    B = x2 - x1
    C = x1 * y2 - x2 * y1

    # Perpendicular distances
    distances = np.abs(A * x + B * y + C) / denom

    # Exclude endpoints
    distances[0] = -np.inf
    distances[-1] = -np.inf

    idx_elbow = int(np.argmax(distances))

    gamma_elbow = gamma_vec[idx_elbow]

    return gamma_elbow

class AnchorPointsEstimation:
    def __init__(self, X, Yt, K, black_box,
                 estimation_method="empirical",
                 gamma = None,
                 calibrate_gamma = False,
                 gamma_vec = None,
                 elbow_detection_method="D2L",
                 drop = 0.01,
                 verbose=False):
        
        self.K = K
        self.method = estimation_method
        self.black_box = copy.deepcopy(black_box)
        self.n_cal = X.shape[0]

        self.calibrate_gamma = calibrate_gamma
        if gamma_vec is not None:
            self.gamma_vec = gamma_vec

        self.elbow_detection_method = elbow_detection_method
        self.drop = drop

        # Estimate the probabilities on the calibration set
        p_hat = black_box.predict_proba(X)
        if not isinstance(p_hat, np.ndarray):
            p_hat = np.asarray(p_hat)

        ## Estimate the contamination process by identifying anchor points
        if self.calibrate_gamma:
            anchor_points_list = {}
            accuracy_tilde_list = []

            for gamma in gamma_vec:
                # Find set of anchor points with threshold gamma%
                anchor_points_list[gamma] = []
                m = max(1, int(np.ceil(gamma * self.n_cal)))
                
                for l in range(self.K):
                    scores_l = p_hat[:, l]
                    top_idx = np.argsort(scores_l)[::-1][:m]
                    # Update the list of anchor points
                    anchor_points_list[gamma].append(top_idx)
                
                accuracy_tilde_gamma = evaluate_accuracy_tilde(Yt, self.K, anchor_points_list[gamma])
                accuracy_tilde_list.append(accuracy_tilde_gamma)

            accuracy_tilde_vec = np.asarray(accuracy_tilde_list)
            self.accuracy_tilde_vec = accuracy_tilde_vec
            if self.elbow_detection_method=="D2L":
                gamma_opt = elbow_gamma_distance_to_line(gamma_vec, accuracy_tilde_vec, smooth_flag=True)
            elif self.elbow_detection_method=="drop":
                gamma_opt = gamma_by_percent_drop(gamma_vec, accuracy_tilde_vec, drop=self.drop)
            self.gamma_opt = gamma_opt
            self.anchor_points = anchor_points_list[gamma_opt].copy()
        else:
            self.gamma_opt = gamma
            # Identify the set of anchor points with threshold gamma%
            anchor_points = []
            m = max(1, int(np.ceil(self.gamma_opt * self.n_cal)))

            for l in range(self.K):
                scores_l = p_hat[:, l]
                top_idx = np.argsort(scores_l)[::-1][:m]
                # Update the list of anchor points
                anchor_points.append(top_idx)

            self.anchor_points = anchor_points

        # Estimate T_hat
        if estimation_method=="empirical":
            T_hat = np.zeros((self.K, self.K), dtype=float)
            for l in range(self.K):
                top_idx = self.anchor_points[l]
                # Use anchor points as true classes and then evaluate empirical frequencies
                if len(top_idx)>0:
                    counts = np.bincount(Yt[top_idx], minlength=self.K)
                    T_hat[:,l] = counts/len(top_idx)
                else:
                    T_hat[:, l] = np.ones(self.K) / self.K
        elif estimation_method=="Patrini":
            T_hat = np.zeros((self.K, self.K), dtype=float)
            for l in range(self.K):
                top_idx = self.anchor_points[l]
                col_T_hat = p_hat[top_idx, :].mean(axis=0)
                T_hat[:, l] = col_T_hat/np.sum(col_T_hat)
        elif estimation_method=="empirical_parametricRR":
            # Evaluate A_tilde in the set of anchor points
            anchors_correct_tilde = 0
            anchors_total = 0
            for l in range(K):
                top_idx = anchor_points[l]
                Yt_top = Yt[top_idx]
                anchors_correct_tilde += np.sum(Yt_top == l)
                anchors_total   += len(top_idx)
            # Use A_tilde as an estimate of epsilon
            accuracy_tilde_gamma = anchors_correct_tilde/anchors_total
            epsilon_hat = K/(K-1)*(1-accuracy_tilde_gamma)
            T_hat = (1-epsilon_hat) * np.identity(K) + epsilon_hat/K * np.ones((K,K))
            
        self.T_hat = T_hat

        if verbose:
            print('Estimation via anchor points completed.')
            sys.stdout.flush()

    def get_estimate(self):
        if self.calibrate_gamma:
            return self.T_hat, self.anchor_points, self.gamma_opt, self.accuracy_tilde_vec
        else:
            return self.T_hat, self.anchor_points, None, None
    
    def get_anchor_dataset(self, X):
        # Construct the anchor dataset
        K = self.K

        p_hat = self.black_box.predict_proba(X)
        if not isinstance(p_hat, np.ndarray):
            p_hat = np.asarray(p_hat)

        # Collect (index, class) pairs
        anchor_pairs = []
        for l in range(K):
            for idx in self.anchor_points[l]:
                anchor_pairs.append((int(idx), l))

        # Resolve possible duplicates:
        # assign each index to the class with maximal p_hat
        anchor_dict = {}
        for idx, l in anchor_pairs:
            if idx not in anchor_dict:
                anchor_dict[idx] = l
            else:
                # choose class with highest predicted probability
                if p_hat[idx, l] > p_hat[idx, anchor_dict[idx]]:
                    anchor_dict[idx] = l

        anchor_idx = np.array(sorted(anchor_dict.keys()), dtype=int)
        Y_anchor = np.array([anchor_dict[i] for i in anchor_idx], dtype=int)
        X_anchor = X[anchor_idx]

        return (X_anchor, Y_anchor)
