import numpy as np
from scipy.special import softmax, log_softmax
from scipy.optimize import minimize
from dataclasses import dataclass
from typing import Optional
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

def contamination_matrix(eps: float, K: int) -> np.ndarray:
    """
    T(eps)_{kl} = (1 - eps) * delta_{kl} + eps / K
    Shape: (K, K).  T[k, l] = P(Y_tilde = k | Y = l).
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
    eps: float,
    K: int,
) -> np.ndarray:
    """
    Compute responsibilities gamma[i, l] = P(Y_i = l | Y~_i, X_i; theta).

    gamma[i, l]  proportional to  T[Y~_i, l] * p_l(X_i; beta)

    Returns gamma of shape (n_c, K).
    """
    # Classifier probabilities: (n_c, K)
    log_p = log_softmax(X_cont @ beta, axis=1)   # log P(Y=l | X_i)

    T = contamination_matrix(eps, K)              # (K, K)

    # log T[Y~_i, l] for each sample i: shape (n_c, K)
    log_T = np.log(T[Y_tilde, :] + 1e-300)       # T[Y~_i, :] selects row

    log_gamma = log_p + log_T                     # (n_c, K)

    # Normalise in log-space for numerical stability
    log_gamma -= log_gamma.max(axis=1, keepdims=True)
    gamma = np.exp(log_gamma)
    gamma /= gamma.sum(axis=1, keepdims=True)
    return gamma                                   # (n_c, K)


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
# M-step: epsilon  (closed form for randomized response)
# ---------------------------------------------------------------------------

def m_step_eps(
    gamma: np.ndarray,
    Y_tilde: np.ndarray,
    K: int,
    eps_bounds: tuple[float, float] = (1e-6, 1.0 - 1e-6),
) -> float:
    """
    Closed-form update for epsilon under the RR contamination model.

    Q(eps) = sum_i sum_l gamma[i,l] * log T(eps)[Y~_i, l]
           = sum_i gamma[i, Y~_i] * log(1 - eps + eps/K)
           + sum_i sum_{l != Y~_i} gamma[i, l] * log(eps/K)

    dQ/d(eps) = 0  =>  closed form below.
    """
    n_c = len(Y_tilde)
    if n_c == 0:
        return eps_bounds[0]

    # Weight on the "diagonal" (Y~ == l) averaged over contaminated samples
    diag_weight = gamma[np.arange(n_c), Y_tilde].mean()

    # Derivation:
    #   eps* = 1 - K/(K-1) * (diag_weight - 1/K)
    eps_new = 1.0 - (K / (K - 1)) * (diag_weight - 1.0 / K)
    return float(np.clip(eps_new, *eps_bounds))


# ---------------------------------------------------------------------------
# Log-likelihood monitor
# ---------------------------------------------------------------------------

def observed_log_likelihood(
    data: Dataset,
    beta: np.ndarray,
    eps: float,
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
        T = contamination_matrix(eps, K)
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
    eps: float
    log_likelihoods: list[float]
    n_iter: int
    converged: bool


def run_em(
    data: Dataset,
    beta_init: Optional[np.ndarray] = None,
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

    eps = float(np.clip(eps_init, 1e-6, 1 - 1e-6))

    # Partition data once
    X_clean = data.X[data.clean_mask]
    Y_clean = data.Y_obs[data.clean_mask]
    X_cont  = data.X[data.cont_mask]
    Y_tilde = data.Y_obs[data.cont_mask]

    log_likelihoods: list[float] = []
    converged = False

    if verbose:
        print(f"{'Iter':>5}  {'Log-lik':>14}  {'Delta':>12}  {'eps':>8}")
        print("-" * 46)

    for t in range(max_iter):
        # ------------------------------------------------------------------
        # E-step: compute responsibilities for contaminated samples
        # ------------------------------------------------------------------
        if len(X_cont):
            gamma = e_step(X_cont, Y_tilde, beta, eps, K)
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
        eps = m_step_eps(gamma, Y_tilde, K)

        # ------------------------------------------------------------------
        # Monitor observed log-likelihood
        # ------------------------------------------------------------------
        ll = observed_log_likelihood(data, beta, eps)
        log_likelihoods.append(ll)

        if verbose:
            delta = ll - log_likelihoods[-2] if t > 0 else float("nan")
            print(f"{t+1:>5}  {ll:>14.4f}  {delta:>12.6f}  {eps:>8.4f}")

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

    return EMResult(
        beta=beta,
        eps=eps,
        log_likelihoods=log_likelihoods,
        n_iter=len(log_likelihoods),
        converged=converged,
    )


# ---------------------------------------------------------------------------
# Prediction helpers
# ---------------------------------------------------------------------------

def predict_proba(X: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Softmax probabilities, shape (n, K)."""
    return softmax(X @ beta, axis=1)


def predict(X: np.ndarray, beta: np.ndarray) -> np.ndarray:
    """Hard class predictions, shape (n,)."""
    return predict_proba(X, beta).argmax(axis=1)