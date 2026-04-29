import numpy as np
from scipy.special import softmax, log_softmax
from scipy.optimize import minimize
from dataclasses import dataclass
from typing import Optional, Literal
import warnings


"""
EM Algorithm for Classification with Contaminated Labels
=========================================================
Model:  Y  ~ Categorical(softmax(X @ beta))   [true label, K classes]
        Y~ ~ T(eps) @ Y                         [contaminated label]
        I_i in {0,1}  (exogenous, known)
          I=1  =>  we observe Y_i   (clean)
          I=0  =>  we observe Y~_i  (contaminated)

Contamination (randomized response):
        T(eps)_{kl} = (1-eps)*I[k==l] + eps/K

Parameters to estimate: beta (logistic weights), eps (contamination rate).
"""

# ---------------------------------------------------------------------------
# Contamination matrix
# ---------------------------------------------------------------------------

def rr_matrix(eps: float, K: int) -> np.ndarray:
    """
    Randomised-response contamination matrix.
    T(eps)[k, l] = (1 - eps) * delta_{kl} + eps / K
    Shape: (K, K).
    """
    return (1.0 - eps) * np.eye(K) + (eps / K) * np.ones((K, K))


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class Dataset:
    X: np.ndarray        # (n, d)  feature matrix
    Y_obs: np.ndarray    # (n,)    observed label (Y if I=1, Y~ if I=0)
    I: np.ndarray        # (n,)    indicator: 1 = clean, 0 = contaminated
    K: int               # number of classes

    @property
    def clean_mask(self) -> np.ndarray:
        return self.I == 1

    @property
    def cont_mask(self) -> np.ndarray:
        return self.I == 0


# ---------------------------------------------------------------------------
# E-step
# ---------------------------------------------------------------------------

def e_step(
    X_cont: np.ndarray,   # (n_c, d)
    Y_tilde: np.ndarray,  # (n_c,)   contaminated labels (0-indexed)
    beta: np.ndarray,     # (d, K)
    T: np.ndarray,        # (K, K)   current contamination matrix
) -> np.ndarray:
    """
    Compute responsibilities:
        gamma[i, l] = P(Y_i = l | Y~_i, X_i; beta, T)
                    proportional to  T[Y~_i, l] * p_l(X_i; beta)

    Returns gamma of shape (n_c, K).
    """
    log_p = log_softmax(X_cont @ beta, axis=1)        # (n_c, K)
    log_T = np.log(T[Y_tilde, :] + 1e-300)            # (n_c, K)  row = Y~_i

    log_gamma = log_p + log_T
    log_gamma -= log_gamma.max(axis=1, keepdims=True)  # numerical stability
    gamma = np.exp(log_gamma)
    gamma /= gamma.sum(axis=1, keepdims=True)
    return gamma                                  # (n_c, K)


# ---------------------------------------------------------------------------
# M-step: beta  (logistic regression on soft labels)
# ---------------------------------------------------------------------------

def _neg_q_beta(
    beta_flat: np.ndarray,
    X_clean: np.ndarray,
    Y_clean: np.ndarray,
    X_cont: np.ndarray,
    gamma: np.ndarray,
    K: int,
    l2_reg: float,
) -> tuple[float, np.ndarray]:
    """Negative Q w.r.t. beta (with gradient) for scipy.optimize."""
    d = X_clean.shape[1] if len(X_clean) else X_cont.shape[1]
    beta = beta_flat.reshape(d, K)

    loss = 0.0
    grad = np.zeros_like(beta)

    # --- Clean observations: standard cross-entropy ---
    if len(X_clean):
        log_p_clean = log_softmax(X_clean @ beta, axis=1)   # (n_cl, K)
        loss -= log_p_clean[np.arange(len(Y_clean)), Y_clean].sum()
        p_clean = np.exp(log_p_clean)
        one_hot = np.zeros_like(p_clean)
        one_hot[np.arange(len(Y_clean)), Y_clean] = 1.0
        grad -= X_clean.T @ (one_hot - p_clean)

    # --- Contaminated observations: soft cross-entropy ---
    if len(X_cont):
        log_p_cont = log_softmax(X_cont @ beta, axis=1)     # (n_c, K)
        loss -= (gamma * log_p_cont).sum()
        p_cont = np.exp(log_p_cont)
        grad -= X_cont.T @ (gamma - p_cont)

    # L2 regularisation
    if l2_reg > 0:
        loss += 0.5 * l2_reg * (beta ** 2).sum()
        grad += l2_reg * beta

    return float(loss), grad.ravel()


def m_step_beta(
    X_clean: np.ndarray,
    Y_clean: np.ndarray,
    X_cont: np.ndarray,
    gamma: np.ndarray,
    K: int,
    beta_init: np.ndarray,
    l2_reg: float = 1e-4,
    max_iter: int = 200,
) -> np.ndarray:
    """
    Maximise Q w.r.t. beta using L-BFGS-B.
    Returns beta of shape (d, K).
    """
    result = minimize(
        fun=_neg_q_beta,
        x0=beta_init.ravel(),
        args=(X_clean, Y_clean, X_cont, gamma, K, l2_reg),
        method="L-BFGS-B",
        jac=True,
        options={"maxiter": max_iter, "ftol": 1e-10, "gtol": 1e-6},
    )
    return result.x.reshape(beta_init.shape)


# ---------------------------------------------------------------------------
# M-step: T  — two variants
# ---------------------------------------------------------------------------

def m_step_T_general(
    gamma: np.ndarray,    # (n_c, K)
    Y_tilde: np.ndarray,  # (n_c,)
    K: int,
    floor: float = 1e-8,
) -> np.ndarray:
    """
    Closed-form M-step for a general row-stochastic T.

    Aggregated responsibility matrix:
        N[k, l] = sum_{i : Y~_i = k}  gamma[i, l]

    Update:
        T[k, l] = N[k, l] / sum_{l'} N[k, l']

    Rows for classes never observed as Y~ are set to uniform as a fallback.
    """
    N = np.zeros((K, K))
    np.add.at(N, Y_tilde, gamma)         # N[Y~_i, :] += gamma[i, :]

    col_sums = N.sum(axis=0, keepdims=True)          # (1, K)
    unobserved = (col_sums.ravel() == 0)
    col_sums = np.where(col_sums == 0, 1.0, col_sums)
    T_new = N / col_sums
    T_new[:, unobserved] = 1.0 / K                  # uniform fallback

    T_new = np.clip(T_new, floor, None)
    T_new /= T_new.sum(axis=0, keepdims=True)        # renormalise columns
    return T_new


def m_step_T_rr(
    gamma: np.ndarray,
    Y_tilde: np.ndarray,
    K: int,
    eps_bounds: tuple[float, float] = (1e-6, 1.0 - 1e-6),
    **kwargs,
) -> np.ndarray:
    """
    Closed-form M-step for the randomised-response model.

        eps* = K/(K-1) * (1 - mean_i gamma[i, Y~_i])

    Returns the full (K, K) T matrix for a uniform interface.
    """
    n_c = len(Y_tilde)
    if n_c == 0:
        return rr_matrix(eps_bounds[0], K)

    diag_weight = gamma[np.arange(n_c), Y_tilde].mean()
    eps = float(np.clip(
        (K / (K - 1)) * (1.0 - diag_weight),
        *eps_bounds,
    ))
    return rr_matrix(eps, K)


# ---------------------------------------------------------------------------
# Log-likelihood monitor
# ---------------------------------------------------------------------------

def observed_log_likelihood(
    data: Dataset,
    beta: np.ndarray,
    T: np.ndarray,
) -> float:
    """
    Observed-data log-likelihood.

    Clean:        sum_{i in O}  log p_{Y_i}(X_i; beta)
    Contaminated: sum_{i in C}  log sum_l T[Y~_i, l] * p_l(X_i; beta)
    """
    K = data.K
    ll = 0.0

    if data.clean_mask.any():
        X_cl = data.X[data.clean_mask]
        Y_cl = data.Y_obs[data.clean_mask]
        log_p = log_softmax(X_cl @ beta, axis=1)
        ll += log_p[np.arange(len(Y_cl)), Y_cl].sum()

    if data.cont_mask.any():
        X_co = data.X[data.cont_mask]
        Y_ti = data.Y_obs[data.cont_mask]
        p = softmax(X_co @ beta, axis=1)                # (n_c, K)
        marginal = (T[Y_ti, :] * p).sum(axis=1)         # P(Y~ | X)
        ll += np.log(marginal + 1e-300).sum()

    return float(ll)


# ---------------------------------------------------------------------------
# Full EM loop
# ---------------------------------------------------------------------------

@dataclass
class EMResult:
    beta: np.ndarray
    T: np.ndarray
    log_likelihoods: list[float]
    n_iter: int
    converged: bool

    @property
    def eps(self) -> float:
        """
        Implied epsilon under the RR parameterisation.
        For a general T this is the average off-diagonal mass (approximate).
        """
        K = self.T.shape[0]
        diag_mean = np.diag(self.T).mean()
        return float(np.clip((K / (K - 1)) * (1.0 - diag_mean), 0.0, 1.0))


def run_em(
    data: Dataset,
    contamination_model_: Literal["uniform", "general"] = "uniform",
    beta_init: Optional[np.ndarray] = None,
    T_init: Optional[np.ndarray] = None,
    eps_init: float = 0.3,
    max_iter: int = 100,
    tol: float = 1e-6,
    l2_reg: float = 1e-4,
    verbose: bool = True,
) -> EMResult:
    """
    Run the EM algorithm.

    Parameters
    ----------
    data      : Dataset object
    beta_init : initial (d, K) weight matrix; random if None
    eps_init  : initial contamination parameter in (0,1)
    max_iter  : maximum number of EM iterations
    tol       : convergence threshold on |Delta log-lik|
    l2_reg    : L2 penalty on beta for M-step
    verbose   : print progress

    Returns
    -------
    EMResult with fitted parameters and diagnostics
    """
    n, d = data.X.shape
    K = data.K

    # Initialise beta
    if beta_init is None:
        rng = np.random.default_rng(42)
        beta = rng.standard_normal((d, K)) * 0.01
    else:
        beta = beta_init.copy()

    # Initialise T
    if T_init is not None:
        T = T_init.copy()
    else:
        T = rr_matrix(eps_init, K)
    
    # Select M-step for T
    if contamination_model_ == "uniform":
        _m_step_T = m_step_T_rr
    elif contamination_model_ == "general":
        _m_step_T = m_step_T_general
    else:
        raise ValueError(f"Unknown contamination model in EM: {contamination_model_!r}")

    # Partition data once
    X_clean = data.X[data.clean_mask]
    Y_clean = data.Y_obs[data.clean_mask]
    X_cont  = data.X[data.cont_mask]
    Y_tilde = data.Y_obs[data.cont_mask]

    log_likelihoods: list[float] = []
    converged = False

    if verbose:
        header = f"{'Iter':>5}  {'Log-lik':>14}  {'Delta':>12}"
        if contamination_model_ == "uniform":
            header += f"  {'eps':>8}"
        print(header)
        print("-" * (len(header) + 2))

    for t in range(max_iter):
        # ------------------------------------------------------------------
        # E-step: compute responsibilities for contaminated samples
        # ------------------------------------------------------------------
        if len(X_cont):
            gamma = e_step(X_cont, Y_tilde, beta, T)
        else:
            gamma = np.empty((0, K))

        # ------------------------------------------------------------------
        # M-step: update beta
        # ------------------------------------------------------------------
        beta = m_step_beta(
            X_clean, Y_clean, X_cont, gamma, K, beta, l2_reg=l2_reg
        )

        # ------------------------------------------------------------------
        # M-step: update eps (closed form)
        # ------------------------------------------------------------------
        T = _m_step_T(gamma, Y_tilde, K)

        # ------------------------------------------------------------------
        # Monitor observed log-likelihood
        # ------------------------------------------------------------------
        ll = observed_log_likelihood(data, beta, T)
        log_likelihoods.append(ll)

        if verbose:
            delta = ll - log_likelihoods[-2] if t > 0 else float("nan")
            row = f"{t+1:>5}  {ll:>14.4f}  {delta:>12.6f}"
            if contamination_model_ == "uniform":
                row += f"  {(K/(K-1))*(1-np.diag(T).mean()):>8.4f}"
            print(row)

        # ------------------------------------------------------------------
        # Convergence check
        # ------------------------------------------------------------------
        if t > 0 and abs(ll - log_likelihoods[-2]) < tol:
            if verbose:
                print(f"\nConverged at iteration {t+1}.")
            converged = True
            break
    else:
        if verbose:
            warnings.warn(f"EM did not converge in {max_iter} iterations.")

    return EMResult(beta=beta, T=T, log_likelihoods=log_likelihoods,
                    n_iter=len(log_likelihoods), converged=converged)


# ---------------------------------------------------------------------------
# Prediction helpers
# ---------------------------------------------------------------------------

def predict_proba(X: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Softmax probabilities, shape (n, K)."""
    return softmax(X @ beta, axis=1)


def predict(X: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Hard class predictions, shape (n,)."""
    return predict_proba(X, beta).argmax(axis=1)