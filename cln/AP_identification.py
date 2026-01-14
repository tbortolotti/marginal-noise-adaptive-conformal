import numpy as np
from sklearn.covariance import LedoitWolf, OAS
from sklearn.mixture import GaussianMixture
import copy

def evaluate_accuracy_tilde(Yt, K, anchor_points_):
    # Evaluate quality of the set of anchor points
    anchors_correct_tilde = 0
    anchors_total = 0

    for l in range(K):
        top_idx = anchor_points_[l]
        Yt_top = Yt[top_idx]
        anchors_correct_tilde += np.sum(Yt_top == l)
        anchors_total += len(top_idx)
    
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

def mahalanobis_dist_(X, mean_, L):
    Xc = X - mean_
    u = np.linalg.solve(L, Xc.T).T
    return np.sqrt(np.sum(u * u, axis=1))

# AP filter based on Mahalanobis depth
def mahalanobis_depth_certifier(X_, Y_, cov_estimator="ledoitwolf", anchor_frac=0.5, anchor_m=None, eps=1e-9):
    
    n_ = X_.shape[0]
    p_ = X_.shape[1]
    classes_ = np.unique(Y_)
    K_ = len(classes_)

    # Compute covariance and cholesky factors
    means_ = np.zeros((K_, p_), dtype=float)
    chol_factors_ = []
    for i, k in enumerate(classes_):
        Xk = X_[Y_ == k]
        means_[i] = Xk.mean(axis=0)

        if cov_estimator.lower() == "empirical":
            cov = np.cov(Xk, rowvar=False, bias=False)
        elif cov_estimator.lower() == "oas":
            cov = OAS().fit(Xk).covariance_
        else:
            cov = LedoitWolf().fit(Xk).covariance_

        cov = cov + eps * np.eye(p_)
        L = np.linalg.cholesky(cov)
        chol_factors_.append(L)

    dists = np.zeros((n_, K_), dtype=float)
    for i in range(K_):
        dists[:, i] = mahalanobis_dist_(X_, means_[i], chol_factors_[i])

    depths = 1.0 / (1.0 + dists)

    top = np.argmax(depths, axis=1)
    top_best = depths[np.arange(n_), top]

    Y_hat = np.full(n_, -1, dtype=int)

    for class_idx, class_label in enumerate(classes_):
        idx_k = np.where(top == class_idx)[0]
        if idx_k.size == 0:
            continue

        # choose how many anchors to keep in this class
        if anchor_m is not None:
            m_ = int(anchor_m)
            m_ = max(0, min(m_, idx_k.size))
        else:
            frac = float(anchor_frac)
            if not (0.0 < frac <= 1.0):
                raise ValueError("anchor_frac must be in (0, 1].")
            m_ = int(np.ceil(frac * idx_k.size))
            m_ = max(1, min(m_, idx_k.size))

        # sort by depth (descending) within this predicted class
        order = np.argsort(-top_best[idx_k])
        keep = idx_k[order[:m_]]

        Y_hat[keep] = class_label

    #Y_hat = np.where((top_best >= self.depth_min) & (gap >= self.depth_gap), self.classes_[top], -1)

    return Y_hat

"""
def multimodal_gmm_certifier(
    X_, Y_,
    n_components=2, # int or dict: {class_label: n_components}
    covariance_type="full", # "full", "diag", "tied", "spherical"
    reg_covar=1e-6,
    anchor_frac=0.5,
    anchor_m=None,
    random_state=0,
    score_mode="class_log_density"  # "class_log_density" or "max_component"
):

    #Multimodal generalization of MahalanobisDepthCertifier using per-class GMMs.

    n_ = X_.shape[0]
    classes_ = np.unique(Y_)
    K_ = len(classes_)

    # Fit one GMM per class
    gmms = {}
    for k in classes_:
        Xk = X_[Y_ == k]
        nk = Xk.shape[0]
        # choose components for class k
        if isinstance(n_components, dict):
            mk = int(n_components[k])
        else:
            mk = int(n_components)

        mk = max(1, min(mk, nk))  # can't exceed number of points

        gmm = GaussianMixture(
            n_components=mk,
            covariance_type=covariance_type,
            reg_covar=reg_covar,
            random_state=random_state
        )
        gmm.fit(Xk)
        gmms[k] = gmm

    # Compute class scores for every point
    # score_samples returns log p(x)
    class_scores = np.zeros((n_, K_), dtype=float)
    for j, k in enumerate(classes_):
        gmm = gmms[k]
        if score_mode == "class_log_density":
            class_scores[:, j] = gmm.score_samples(X_)   # log p_k(x)
        elif score_mode == "max_component":
            # responsibilities proportional to pi_m N_m(x); use component-wise log prob
            # sklearn doesn't expose per-component log N directly, so we compute log prob of each component:
            # Using internal helper is not public; simplest robust approach:
            # approximate max component score via logsumexp trick + responsibilities is more complex
            # -> alternative: use predict_proba responsibilities and back out per-component unnormalized log probs
            resp = gmm.predict_proba(X_)  # shape (n, mk)
            # take max responsibility as a proxy for being near one mode (works well in practice)
            class_scores[:, j] = resp.max(axis=1)
        else:
            raise ValueError("score_mode must be 'class_log_density' or 'max_component'.")

    # Predicted class = argmax class score
    top = np.argmax(class_scores, axis=1)
    top_best = class_scores[np.arange(n_), top]

    Y_hat = np.full(n_, -1, dtype=int)

    for class_idx, class_label in enumerate(classes_):
        idx_k = np.where(top == class_idx)[0]
        if idx_k.size == 0:
            continue

        # choose how many anchors to keep in this predicted class
        if anchor_m is not None:
            m_ = int(anchor_m)
            m_ = max(0, min(m_, idx_k.size))
        else:
            frac = float(anchor_frac)
            if not (0.0 < frac <= 1.0):
                raise ValueError("anchor_frac must be in (0, 1].")
            m_ = int(np.ceil(frac * idx_k.size))
            m_ = max(1, min(m_, idx_k.size))

        order = np.argsort(-top_best[idx_k])  # descending score
        keep = idx_k[order[:m_]]
        Y_hat[keep] = class_label

    return Y_hat
"""



class AnchorPointsIdentification:
    def __init__(self, X, Yt, K, black_box,
                 gamma = None,
                 calibrate_gamma = False,
                 gamma_vec = None,
                 elbow_detection_method="D2L",
                 drop = 0.01,
                 ap_filter=False,
                 filter_method="mahalanobisdepth"):
        
        self.K = K
        self.classes_ = np.unique(Yt)
        self.black_box = copy.deepcopy(black_box)
        self.n_cal = X.shape[0]

        self.calibrate_gamma = calibrate_gamma
        if gamma_vec is not None:
            self.gamma_vec = gamma_vec

        self.elbow_detection_method = elbow_detection_method
        self.drop = drop

        self.ap_filter = ap_filter
        self.filter_method = filter_method

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
        
        if ap_filter:
            # Operate a filtering of the anchor points set in order to improve accuracy
            if filter_method=="mahalanobisdepth":
                X_anchor, Yt_anchor, anchor_idx = self.get_anchor_dataset(X)
                Yt_anchor_filtered = mahalanobis_depth_certifier(X_anchor, Yt_anchor, cov_estimator="ledoitwolf", anchor_m=20)

                keep_mask = (Yt_anchor_filtered != -1)
                kept_indices = set(anchor_idx[keep_mask])

                for l in range(self.K):
                    self.anchor_points[l] = [
                        idx for idx in self.anchor_points[l] if idx in kept_indices
                    ]

    def get_ap_(self):
        if self.calibrate_gamma:
            return self.anchor_points, self.gamma_opt, self.accuracy_tilde_vec
        else:
            return self.anchor_points, None, None
    
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

        return X_anchor, Y_anchor, anchor_idx

