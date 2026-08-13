"""
Conformal prediction under label contamination, adapted from Clarkson,
Xu, Yang, Zeng (2024), "Conformal Prediction under Data Contamination"
(https://arxiv.org/abs/2407.07700), Section 4.4 / Theorem 5.1.

Ported from the authors' repository (cp_under_data_contamination-master:
robust_cp.py, aps.py, run_experiment.py) into this repository's classification
scaffolding, so that it can be dropped into the same `methods` dictionary
pattern used in experiments/exp_1.py alongside
cln.classification.MarginalLabelNoiseConformal and
cln.classification_label_conditional.LabelNoiseConformal.

No function already available in this repository is duplicated:
 - the noise model is parametrized with the T/M/rho/rho_tilde quantities
   already produced by cln.contamination and cln.data, instead of the
   original repository's separate P/p/p_tilde convention. The two are related
   by P = M.T (see ClarksonConformal.__init__ for the derivation), so no
   information needs to be recomputed from scratch.

Only robust_cp.py's bias-correction machinery (ScoreECDF, g_hat, compute_b,
compute_C, compute_adjusted_i) is genuinely new relative to this repository
and is ported here (Section 4.4 / Theorem 5.1 of the paper).
"""

import os
import sys
import copy

import numpy as np
from sklearn.model_selection import train_test_split
from statsmodels.distributions.empirical_distribution import ECDF

# Make `arc` (and, for the __main__ demo below, `cln`) importable regardless
# of the caller's working directory, mirroring the sys.path setup used in
# experiments/exp_1.py.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_THIRD_PARTY_DIR = os.path.join(_THIS_DIR, "..")
_REPO_ROOT_DIR = os.path.join(_THIS_DIR, "..", "..")
for _path in (_THIRD_PARTY_DIR, _REPO_ROOT_DIR):
    if _path not in sys.path:
        sys.path.append(_path)

from arc.classification import ProbabilityAccumulator as ProbAccum


class ScoreECDF:
    """
    Empirical CDFs of the conformity scores, computed conditionally on each
    possible value of the *contaminated* label, for every candidate label
    used as a placeholder (Section 4.4 of the paper).

    Adapted from robust_cp.py's ScoreECDF: uses this repository's
    third_party.arc.classification.ProbabilityAccumulator (functionally
    identical to the original repository's local aps.ProbabilityAccumulator,
    both being copies of msesia/arc), and, like
    cln.classification.MarginalLabelNoiseConformal.compute_F_hat_scores,
    falls back to the unconditional (global) score distribution for any
    contaminated-label class with zero calibration samples instead of
    leaving a `None` ECDF.
    """

    def __init__(self, Y: np.ndarray, K: int, grey_box: ProbAccum, random_noise: np.ndarray):
        self.K = K
        self.ecdfs = [[None] * K for _ in range(K)]
        self.scores = [[None] * K for _ in range(K)]

        # Conformity scores obtained by pretending the true label is k, for
        # every calibration point (reused across the conditioning classes j).
        global_scores = {}
        for k in range(K):
            Y_k = k * np.ones_like(Y)
            # NOTE: the APS-style score is 1 - (calibrated p-value), matching
            # the "1.0 - ..." convention in robust_cp.py.
            global_scores[k] = 1.0 - grey_box.calibrate_scores(Y_k, epsilon=random_noise)

        for j in range(K):
            idx = np.where(Y == j)[0]
            for i in range(K):
                if len(idx) > 0:
                    self.scores[i][j] = global_scores[i][idx]
                else:
                    self.scores[i][j] = global_scores[i]
                self.ecdfs[i][j] = ECDF(self.scores[i][j])

    def compute(self, q: np.ndarray, i: int, j: int) -> np.ndarray:
        return self.ecdfs[i][j](q)


def g_hat(q: np.ndarray, fhat: ScoreECDF, P_inv: np.ndarray, rho: np.ndarray, rho_tilde: np.ndarray) -> np.ndarray:
    """
    Plug-in estimator of the bias g(q), Equation (31) of the paper.

    P_inv is the inverse of P[j, i] = P(Y=j | Ytilde=i), which, with this
    repository's mixing matrix M (M[k, l] = P(Y=l | Ytilde=k), see
    cln.contamination.convert_T_to_M), equals M.T.
    """
    K = P_inv.shape[0]
    ret = np.zeros_like(q, dtype=float)
    for i in range(K):
        for j in range(K):
            ret += rho[i] * P_inv[j, i] * fhat.compute(q, i, j)
        ret -= rho_tilde[i] * fhat.compute(q, i, i)
    return ret


def compute_b(n: int, rho_tilde: np.ndarray) -> np.ndarray:
    """Coefficients b(n, i) of Theorem 5.1 of the paper."""
    return np.sqrt(np.pi / (n * rho_tilde)) + (1 - rho_tilde) ** n


def compute_C(n: int, P_inv: np.ndarray, rho: np.ndarray, rho_tilde: np.ndarray) -> float:
    """Upper bound on E_q[g_hat(q) - g(q)], Theorem 5.1 of the paper."""
    K = P_inv.shape[0]
    b = compute_b(n, rho_tilde)
    ret = 0.0
    for i in range(K):
        ret += np.abs(P_inv[i, i] * rho[i] - rho_tilde[i]) * b[i]
        for j in range(K):
            if j != i:
                ret += np.abs(P_inv[j, i] * rho[i]) * b[j]
    return ret


def compute_adjusted_i(alpha: float, sorted_scores: np.ndarray, fhat: ScoreECDF,
                        P_inv: np.ndarray, rho: np.ndarray, rho_tilde: np.ndarray,
                        verbose: bool = False) -> int:
    """Adjusted index i into `sorted_scores`, Equation (32) of the paper."""
    n = sorted_scores.shape[0]
    C = compute_C(n, P_inv, rho, rho_tilde)
    ghats = g_hat(sorted_scores, fhat, P_inv, rho, rho_tilde)
    ivals = (1 - alpha) - ghats + C - (np.arange(1, n + 1) / (n + 1))
    candidates = np.where(ivals <= 0)[0]
    if candidates.size > 0:
        i = np.min(candidates)
    else:
        # Fall back to the standard (unadjusted) split-conformal quantile,
        # as in the original repository. Clipped to a valid 0-based index
        # (the original code does not clip, which can go out of bounds).
        i = min(int(np.ceil((n + 1) * (1 - alpha))), n - 1)
        if verbose:
            print("Warning: no valid Clarkson-adjusted index found, "
                  "falling back to standard split-conformal quantile.")
    return i


class ClarksonConformal:
    """
    Marginal conformal prediction under label contamination (Clarkson et al.,
    2024). Same constructor/predict interface as
    cln.classification.MarginalLabelNoiseConformal, so it can be dropped into
    the `methods` dictionary of experiments/exp_1.py, e.g.:

        "Clarkson": lambda: ClarksonConformal(X, Yt, black_box_pt, K, alpha,
                                              n_cal=n_cal, M=M, rho_tilde=rho_tilde_hat,
                                              allow_empty=allow_empty, pre_trained=True,
                                              random_state=random_state),
    """

    def __init__(self, X, Y, black_box, K, alpha, n_cal=0.5, M=None, rho_tilde=None,
                 allow_empty=False, verbose=False, pre_trained=False, random_state=2023):

        self.K = K
        self.allow_empty = allow_empty
        self.black_box = copy.deepcopy(black_box)

        if M is None or rho_tilde is None:
            raise ValueError("Both the mixing matrix M and rho_tilde must be provided.")
        self.M = M
        self.rho_tilde = rho_tilde
        # P[j, i] = P(Y=j | Ytilde=i) in the paper's notation equals M.T here.
        self.P_inv = np.linalg.inv(M.T)
        # Marginal label proportions implied by M and rho_tilde, as in
        # cln.classification.MarginalLabelNoiseConformal.compute_delta_const_marginal.
        self.rho = np.dot(M.T, rho_tilde)

        # Split data into training/calibration sets
        if n_cal >= 0:
            X_train, X_calib, Y_train, Y_calib = train_test_split(X, Y, test_size=n_cal, random_state=random_state)
        else:
            X_calib, Y_calib = X, Y
        n_cal = X_calib.shape[0]

        # Fit model (if not pre-trained)
        if pre_trained:
            if verbose:
                print('Skipping training.')
                sys.stdout.flush()
        else:
            if verbose:
                print('Fitting classifier on {:d} training samples with {:d} features...'.format(
                    X_train.shape[0], X_train.shape[1]))
                sys.stdout.flush()
            self.black_box.fit(X_train, Y_train)
            if verbose:
                print('Training completed.')
                sys.stdout.flush()

        if verbose:
            print('Evaluating conformity scores on {:d} calibration samples...'.format(X_calib.shape[0]))
            sys.stdout.flush()

        # Evaluate conformity scores on calibration data
        p_hat_calib = self.black_box.predict_proba(X_calib)
        grey_box = ProbAccum(p_hat_calib)
        rng = np.random.default_rng(random_state)
        random_noise = rng.uniform(low=0.0, high=1.0, size=n_cal)

        if verbose:
            print('Calibrating conformity scores for {:d} classes...'.format(K))
            sys.stdout.flush()

        fhat = ScoreECDF(Y_calib, K, grey_box, random_noise)
        scores_sorted = np.sort(np.concatenate([fhat.scores[k][k] for k in range(K)]))

        i_adjusted = compute_adjusted_i(alpha, scores_sorted, fhat, self.P_inv, self.rho, self.rho_tilde,
                                        verbose=verbose)
        qhat = scores_sorted[i_adjusted]

        # A single marginal threshold, stored as a length-K array for
        # interface-compatibility with MarginalLabelNoiseConformal.predict().
        self.tau_hat = qhat * np.ones((K,))

        if verbose:
            print('Calibration of conformity scores completed.')
            sys.stdout.flush()

    def predict(self, X, random_state=2023):
        n = X.shape[0]
        rng = np.random.default_rng(random_state)
        random_noise = rng.uniform(low=0.0, high=1.0, size=n)
        p_hat = self.black_box.predict_proba(X)
        grey_box = ProbAccum(p_hat)
        S_hat = [np.array([]).astype(int) for _ in range(n)]
        for k in range(self.K):
            alpha_k = 1.0 - self.tau_hat[k]
            S_k = grey_box.predict_sets(alpha_k, epsilon=random_noise, allow_empty=self.allow_empty)
            for i in range(n):
                if k in S_k[i]:
                    S_hat[i] = np.append(S_hat[i], k)
        return S_hat
