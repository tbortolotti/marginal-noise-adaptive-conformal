import numpy as np
#from sklearn.covariance import LedoitWolf, OAS
from sklearn.mixture import GaussianMixture
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import EllipticEnvelope
from scipy.stats import norm

import copy
import sys

#_________________________________________________________________________________________
# METHODS FOR AUTOMATIC SELECTION OF THE THRESHOLD IN ANCHOR POINTS IDENTIFICATION
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


## _________________________________________________________________________________________
# ANCHOR POINTS IDENTIFICATION VIA INLIER DETECTION
def _fit_detector(X_train, method,
                  n_estimators, n_neighbors, random_state=None, max_features=None):
    """
    Returns anomaly scores s for X_test, where larger = more anomalous.
    """
    if method == "elliptic_envelope":
        clf = EllipticEnvelope()
        clf.fit(X_train)

    elif method == "isolation_forest":
        if max_features is None:
            max_features = X_train.shape[1]
        clf = IsolationForest(n_estimators=n_estimators, max_features=max_features,random_state=random_state)
        clf.fit(X_train)

    elif method == "lof":
        clf = LocalOutlierFactor(n_neighbors=n_neighbors, novelty=True)
        clf.fit(X_train)

    else:
        raise ValueError(f"Unknown outlier_detection_method='{method}'")

    return clf

def _select_inliers_from_scores(scores, selection="fixed", inlier_frac=0.8, gmm_random_state=None):
    """
    Convert anomaly scores to an inlier mask.

    scores: array, larger = more anomalous
    selection:
      - "fixed": keep the lowest `inlier_frac` fraction
      - "accuracy": select the optimal threshold at the elbow of accuracy
      - "gmm": fit 2-component GMM; keep points assigned to the low-mean component
    """

    if selection == "fixed":
        if not (0.0 < inlier_frac <= 1.0):
            raise ValueError("inlier_frac must be in (0, 1].")
        thr = np.quantile(scores, inlier_frac)
        inlier_mask = scores <= thr
        return inlier_mask, thr, inlier_mask.mean()

    if selection == "gmm":
        # Data-driven split: two groups with different score distributions
        # Inliers are expected to have LOWER anomaly scores.
        gmm = GaussianMixture(n_components=2, random_state=gmm_random_state)
        gmm.fit(scores.reshape(-1, 1))

        means = gmm.means_.ravel()
        inlier_comp = np.argmin(means)

        post = gmm.predict_proba(scores.reshape(-1, 1))[:, inlier_comp]
        inlier_mask = post >= 0.5

        # A “threshold” in score space is not unique under GMM, but we can report the
        # max score among selected inliers as an operational cutoff.
        thr = np.max(scores[inlier_mask]) if np.any(inlier_mask) else np.min(scores)
        return inlier_mask, thr, inlier_mask.mean()

    raise ValueError(f"Unknown selection='{selection}'")


def outlier_detection_(X1_, Y1_, X2_, Y2_, K_,
                       outlier_detection_method="isolation_forest",
                       n_estimators=100, n_neighbors=5, max_features=None,
                       selection="fixed",
                       threshold_vec=None,
                       inlier_frac_vec=None,
                       random_state=None):
    """
    Fits a separate detector for each class k using (X1_k, Y1==k),
    scores the points in X2_k, and flags as outliers according to `selection`.

    Output:
      Y2_hat: copy of Y2_ where outliers are set to -1.
    """
    n2_ = X2_.shape[0]
    Y2_hat = -np.ones(n2_, dtype=int)
    scores = np.full((n2_, K_), np.nan, dtype=float)

    if threshold_vec is None:
        threshold_vec = np.ones(K_) / K_

    for k in range(K_):
        idx_train = (Y1_ == k)
        X1_k = X1_[idx_train]

        # Score X2_k (higher score = more anomalous)
        clf_k = _fit_detector(X_train=X1_k,
                            method=outlier_detection_method,
                            n_estimators=n_estimators, n_neighbors=n_neighbors, max_features=max_features,
                            random_state=random_state)
        sk = -clf_k.score_samples(X2_)
        sk = np.asarray(sk).reshape(-1)
        if sk.shape[0] != n2_:
            raise ValueError(f"Detector {k} score_samples returned shape {sk.shape}, expected ({n2_},)")
        scores[:, k] = sk

    # ---------- Anchor-point selection starts here ----------
    # Require that a point has finite scores from all K detectors.
    # If some detectors are missing (NaN columns), we can only decide among available ones.
    finite = np.isfinite(scores)

    # If a whole row has <2 finite entries, it cannot be "low for one and high for all others".
    row_finite_counts = finite.sum(axis=1)
    valid_row = row_finite_counts == K_

    if selection=="accuracy":

        size_vec = np.zeros(len(inlier_frac_vec))
        accuracy_tilde_vec = np.zeros(len(inlier_frac_vec))
        accuracy_tilde_vec_lower = np.zeros(len(inlier_frac_vec))
        accuracy_tilde_vec_upper = np.zeros(len(inlier_frac_vec))

        for i, inlier_frac_ in enumerate(inlier_frac_vec):
            Y2_hat_ = -np.ones(n2_, dtype=int)
            threshold_vec_i = inlier_frac_ * np.ones(K_)
            inlier_thr = np.full(K_, np.nan)
            for k in range(K_):
                col = scores[:, k]
                col = col[np.isfinite(col)]
                if col.size == 0:
                    continue
                inlier_thr[k] = np.quantile(col, threshold_vec_i[k])

            # Condition: An observation is anchor if it is inlier for exactly one class
            is_inlier = scores <= inlier_thr[None, :]      # (n2_, K_)
            unique_inlier = (is_inlier.sum(axis=1) == 1)
            anchor = valid_row & unique_inlier
            Y2_hat_[anchor] = np.argmin(scores[anchor, :], axis=1)

            accuracy_tilde_vec[i] = np.sum(Y2_hat_==Y2_)/np.sum(Y2_hat_ != -1)
            size_vec[i] = np.sum(Y2_hat_!=-1)

            # Build confidence interval for proportion under Gaussian approximation
            alpha = 0.05
            z = norm.ppf(1 - alpha / 2)
            ci_method="wald"
            if ci_method=="wald":
                se = np.sqrt(accuracy_tilde_vec[i] * (1 - accuracy_tilde_vec[i]) / size_vec[i])
                accuracy_tilde_vec_lower[i] = max(0.0, accuracy_tilde_vec[i] - z * se)
                accuracy_tilde_vec_upper[i] = min(1.0, accuracy_tilde_vec[i] + z * se)
            elif ci_method=="wilson":
                denom = 1 + z*z/size_vec[i]
                center = (accuracy_tilde_vec[i] + z*z/(2*size_vec[i])) / denom
                half = (z / denom) * np.sqrt(accuracy_tilde_vec[i]*(1-accuracy_tilde_vec[i])/size_vec[i] + z*z/(4*size_vec[i]*size_vec[i]))
                accuracy_tilde_vec_lower[i] = max(0.0, center - half)
                accuracy_tilde_vec_upper[i] = min(1.0, center + half)


        #idx = len(accuracy_tilde_vec) - 1 - np.argmax(accuracy_tilde_vec[::-1])
        idx = np.argmax(accuracy_tilde_vec_lower)
        inlier_frac_opt = inlier_frac_vec[idx]
        threshold_vec = inlier_frac_opt * np.ones(K_)

    inlier_thr = np.full(K_, np.nan)
    for k in range(K_):
        col = scores[:, k]
        col = col[np.isfinite(col)]
        if col.size == 0:
            continue
        inlier_thr[k] = np.quantile(col, threshold_vec[k])

        # Condition: An observation is anchor if it is inlier for exactly one class
        is_inlier = scores <= inlier_thr[None, :]      # (n2_, K_)
        unique_inlier = (is_inlier.sum(axis=1) == 1)
        anchor = valid_row & unique_inlier

    # NOTE: the following code is needed if one wants to add a margin
    # Winner (smallest score) and runner-up for margin test
    #best_k = np.argmin(scores, axis=1)
    #best_s = scores[np.arange(n2_), best_k]
    # second smallest score per row
    #second_s = np.partition(scores, 1, axis=1)[:, 1]
    #gaps = second_s - best_s
    #margin = np.quantile(gaps[np.isfinite(gaps)], 0.95)
    #separation_ok = gaps >= margin
    #anchor = valid_row & separation_ok & unique_inlier

    Y2_hat[anchor] = np.argmin(scores[anchor, :], axis=1)

    return Y2_hat, scores

## _________________________________________________________________________________________
# ANCHOR POINTS INDENTIFICATION CLASS
class AnchorPointsIdentification:
    def __init__(self, X1, Yt1, X2, Yt2, K,
                 use_classifier = False,
                 black_box = None,
                 gamma = None,
                 calibrate_gamma = False,
                 gamma_vec = None,
                 outlier_detection=False,
                 outlier_detection_method="isolation_forest",
                 n_neighbors=100,
                 selection="fixed",
                 threshold_vec=None,
                 inlier_frac_vec=None):
        
        # This function produces Y_anchor.
        # Y_anchor is a vector of length equal to the length of Y.
        # It has -1 in correspondence of non-anchor observations, and has value k if the observation is anchor for class k.

        # Strategies for anchor points identification can either be used separately or combined together by conveniently flagging options

        # This class identifies anchor points
        self.K = K
        self.classes_ = np.unique(Yt1)
        self.n2 = X2.shape[0]

        # Set flags that identify AP identification method
        self.use_classifier = use_classifier

        self.calibrate_gamma = calibrate_gamma
        if gamma_vec is not None:
            self.gamma_vec = gamma_vec
        else:
            self.gamma_vec = np.arange(0.1, 0.89 + 1e-12, 0.05)

        self.outlier_detection = outlier_detection
        self.outlier_detection_method = outlier_detection_method
        self.n_neighbors=n_neighbors
        self.selection = selection

        if threshold_vec is None:
            self.threshold_vec = 1/K * np.ones(K)
        else:
            self.threshold_vec = threshold_vec

        if inlier_frac_vec is None:
            rho_k = self.threshold_vec[0]
            min_log = 3
            log_part = np.logspace(-min_log, -1.0, min_log)
            if rho_k > 0.1:
                lin_part = np.arange(0.1 + 0.05, rho_k + 1e-12, 0.05)
                self.inlier_frac_vec = np.concatenate([log_part, lin_part])
            else:
                self.inlier_frac_vec = log_part[log_part <= rho_k]
        else:
            self.inlier_frac_vec=inlier_frac_vec

        if self.use_classifier:
            # Fit the point predictor on the training set
            black_box_pt = copy.deepcopy(black_box)
            black_box_pt.fit(X1, Yt1)

            # Calculate the probabilities on the set for which you want to identify anchor points
            p_hat = black_box_pt.predict_proba(X2)
            if not isinstance(p_hat, np.ndarray):
                p_hat = np.asarray(p_hat)
            self.p_hat = p_hat

            # Estimate the contamination process by identifying anchor points
            if self.calibrate_gamma:
                alpha = 0.05
                size_vec = np.zeros(len(self.gamma_vec))
                accuracy_tilde_vec = np.zeros(len(self.gamma_vec))
                accuracy_tilde_vec_upper = np.zeros(len(self.gamma_vec))
                accuracy_tilde_vec_lower = np.zeros(len(self.gamma_vec))

                for i, gamma in enumerate(self.gamma_vec):
                    Y_anchor_gamma = self.new_identify_anchor_points(gamma)
                    accuracy_tilde_vec[i] = np.sum(Y_anchor_gamma==Yt2)/np.sum(Y_anchor_gamma!=-1)
                    size_vec[i] = np.sum(Y_anchor_gamma!=-1)

                    # Build confidence interval for proportion under Gaussian approximation
                    z = norm.ppf(1 - alpha / 2)
                    se = np.sqrt(accuracy_tilde_vec[i] * (1 - accuracy_tilde_vec[i]) / size_vec[i])
                    accuracy_tilde_vec_lower[i] = max(0.0, accuracy_tilde_vec[i] - z * se)
                    accuracy_tilde_vec_upper[i] = min(1.0, accuracy_tilde_vec[i] + z * se)

                self.accuracy_tilde_vec = accuracy_tilde_vec

                idx = np.argmax(accuracy_tilde_vec_lower)
                gamma_opt = self.gamma_vec[idx]
                self.gamma_opt = gamma_opt
            else:
                self.gamma_opt = gamma
            
            self.Y_anchor = self.new_identify_anchor_points(self.gamma_opt)
            self.scores = p_hat
        
        if self.outlier_detection:
            if self.use_classifier:
                Yt2_ = self.Y_anchor
            else:
                Yt2_ = Yt2.copy()
            # Operate class-wise outlier detection
            # Parameters
            # n_estimators=100, n_neighbors=20, max_features=None, random_state=None,
            # selection="fixed", inlier_frac=0.8
            Yt2_inliers, scores = outlier_detection_(X1, Yt1, X2, Yt2_, self.K,
                                                    outlier_detection_method=outlier_detection_method,
                                                    n_neighbors=self.n_neighbors,
                                                     selection=self.selection,
                                                     threshold_vec=self.threshold_vec,
                                                     inlier_frac_vec=self.inlier_frac_vec)
            self.Y_anchor = Yt2_inliers
            self.scores = scores

    def identify_anchor_points(self, gamma):
        # Identify the set of anchor points with threshold gamma%
        anchor_points = []
        m = max(1, int(np.ceil(gamma * self.n2)))

        for l in range(self.K):
            scores_l = self.p_hat[:, l]
            top_idx = np.argsort(scores_l)[::-1][:m]
            anchor_points.append(top_idx)

        # For every observation, Y_anchor_ stores the class for which the observation is anchor point. Y_anchor_[idx]=-1 if idx is never an anchor point.
        # In doing so, possibe duplicates are resolved:
        # if an observation is anchor point for more than one class, we keep it as anchor point for the class for which it has the highest predicted probability
        Y_anchor_ = -np.ones(self.n2, dtype=np.int32)
        best_prob = -np.inf * np.ones(self.n2)

        for l in range(self.K):
            idx = anchor_points[l]
            probs = self.p_hat[idx, l]

            mask = probs > best_prob[idx]
            Y_anchor_[idx[mask]] = l
            best_prob[idx[mask]] = probs[mask]

        return Y_anchor_
    
    def new_identify_anchor_points(self, threshold):
        # Identify the set of anchor points with threshold
        # score[i,k] = p_hat[i,k] - sum_{l!=k} p_hat[i,l] = 2*p_hat[i,k] - sum_j p[i,j]
        row_sum = self.p_hat.sum(axis=1, keepdims=True)
        scores = 2.0 * self.p_hat - row_sum

        # best class per observation
        best_k = np.argmax(scores, axis=1)
        best_score = scores[np.arange(self.p_hat.shape[0]), best_k]

        Y_anchor_ = -np.ones(self.n2, dtype=np.int32)
        Y_anchor_[best_score > threshold] = best_k[best_score > threshold]

        return Y_anchor_
    

    def get_anchor_points(self):
        if self.calibrate_gamma:
            return self.Y_anchor, self.gamma_opt, self.accuracy_tilde_vec, self.scores
        else:
            return self.Y_anchor, None, None, self.scores

"""
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
"""
