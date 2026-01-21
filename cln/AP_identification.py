import numpy as np
#from sklearn.covariance import LedoitWolf, OAS
#from sklearn.mixture import GaussianMixture
from sklearn.ensemble import IsolationForest
from sklearn.covariance import EllipticEnvelope, LedoitWolf, OAS
import copy
import sys

def gamma_by_percent_drop(gamma_vec, accuracy_vec, drop=0.03):
    acc_max = np.max(accuracy_vec)
    ok = np.where(accuracy_vec >= (1.0-drop)*acc_max)[0]
    if ok.size == 0:
        idx = int(np.argmax(accuracy_vec))
    else:
        idx = int(ok[-1])
    return gamma_vec[idx]

def elbow_gamma_distance_to_line(gamma_vec, accuracy_vec, smooth_flag=False):
    # Distance-to-line elbow finder in log(gamma) space.

    # Convert to numeric arrays
    gamma_vec = np.asarray(gamma_vec, dtype=float).ravel()
    y = np.asarray(accuracy_vec, dtype=float).ravel()

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

def elliptic_envelope_certifier(X_, Y_, K_):
    n_ = X_.shape[0]

    inliers = -np.ones(n_, dtype=np.int32)
    Y_hat = Y_.copy()

    for k in range(K_):
        idxs = (Y_==k)
        X_k = X_[idxs,]
        clf = EllipticEnvelope(random_state=0)
        clf.fit(X_k)
        inliers[idxs] = clf.predict(X_[idxs])

    Y_hat[inliers == -1] = -1
    return Y_hat

def isolation_forest_certifier(X_, Y_, K_, n_estimators=10, warm_start=True):
    n_ = X_.shape[0]

    inliers = -np.ones(n_, dtype=np.int32)
    Y_hat = Y_.copy()

    for k in range(K_):
        idxs = (Y_==k)
        X_k = X_[idxs,]
        clf = IsolationForest(n_estimators=n_estimators, warm_start=warm_start)
        clf.fit(X_k)
        inliers[idxs] = clf.predict(X_[idxs])

    Y_hat[inliers == -1] = -1
    return Y_hat

class AnchorPointsIdentification:
    def __init__(self, X, Yt, K, black_box,
                 gamma = None,
                 calibrate_gamma = False,
                 gamma_vec = None,
                 elbow_detection_method="D2L",
                 min_flag = False,
                 drop = 0.01,
                 ap_filter=False,
                 filter_method="isolation_forest"):
        
        # This function produces Y_anchor.
        # Y_anchor is a vector of length equal to the length of Y.
        # It has -1 in correspondence of non-anchor observations, and has value k if the observation is anchor for class k.
        
        # This class identifies anchor points
        self.K = K
        self.classes_ = np.unique(Yt)
        self.black_box = copy.deepcopy(black_box)
        self.n = X.shape[0]

        self.calibrate_gamma = calibrate_gamma
        if gamma_vec is not None:
            self.gamma_vec = gamma_vec

        self.elbow_detection_method = elbow_detection_method
        self.drop = drop
        self.min_flag = min_flag

        self.ap_filter = ap_filter
        self.filter_method = filter_method

        # Estimate the probabilities on the calibration set
        p_hat = black_box.predict_proba(X)
        if not isinstance(p_hat, np.ndarray):
            p_hat = np.asarray(p_hat)
        self.p_hat = p_hat

        ## Estimate the contamination process by identifying anchor points
        if self.calibrate_gamma:
            accuracy_tilde_vec = np.zeros(len(gamma_vec))
            for i, gamma in enumerate(gamma_vec):
                Y_anchor_gamma = self.identify_anchor_points(gamma)
                accuracy_tilde_vec[i] = np.sum(Y_anchor_gamma==Yt)/np.sum(Y_anchor_gamma!=-1)
            self.accuracy_tilde_vec = accuracy_tilde_vec

            if self.elbow_detection_method=="D2L":
                gamma_opt = elbow_gamma_distance_to_line(gamma_vec, accuracy_tilde_vec, smooth_flag=True)
            elif self.elbow_detection_method=="drop":
                gamma_opt = gamma_by_percent_drop(gamma_vec, accuracy_tilde_vec, drop=self.drop)
            self.gamma_opt = gamma_opt
        else:
            gamma_opt = gamma

        if self.min_flag:
            self.gamma_opt = min(gamma_opt,50*self.K/self.n)
        else:
            self.gamma_opt = gamma_opt
        
        Y_anchor = self.identify_anchor_points(self.gamma_opt)
        
        if ap_filter:
            # Operate a filtering of the anchor points set in order to improve accuracy
            if filter_method=="elliptic_envelope":
                Y_anchor_filtered = elliptic_envelope_certifier(X, Y_anchor, self.K)
            elif filter_method=="isolation_forest":
                Y_anchor_filtered = isolation_forest_certifier(X, Y_anchor, self.K, n_estimators=10, warm_start=True)
            else:
                print("Unknown outlier detection method!")
                sys.stdout.flush()
                exit(-1)
            self.Y_anchor = Y_anchor_filtered
        else:
            self.Y_anchor = Y_anchor

    def identify_anchor_points(self, gamma):
        # Identify the set of anchor points with threshold gamma%
        anchor_points = []
        m = max(1, int(np.ceil(gamma * self.n)))

        for l in range(self.K):
            scores_l = self.p_hat[:, l]
            top_idx = np.argsort(scores_l)[::-1][:m]
            anchor_points.append(top_idx)

        # For every observation, Y_anchor_ stores the class for which the observation is anchor point. Y_anchor_[idx]=-1 if idx is never an anchor point.
        # In doing so, possibe duplicates are resolved:
        # if an observation is anchor point for more than one class, we keep it as anchor point for the class for which it has the highest predicted probability
        Y_anchor_ = -np.ones(self.n, dtype=np.int32)
        best_prob = -np.inf * np.ones(self.n)

        for l in range(self.K):
            idx = anchor_points[l]
            probs = self.p_hat[idx, l]

            mask = probs > best_prob[idx]
            Y_anchor_[idx[mask]] = l
            best_prob[idx[mask]] = probs[mask]

        return Y_anchor_

    def get_anchor_points(self):
        if self.calibrate_gamma:
            return self.Y_anchor, self.gamma_opt, self.accuracy_tilde_vec
        else:
            return self.Y_anchor, None, None
