import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from dataclasses import dataclass, field
from typing import Optional, Literal
import warnings

"""
EM Algorithm for Classification with Contaminated Labels
— Classifier: arbitrary nn.Module (default: MLP) —
=========================================================
This is a drop-in replacement for T_estimation_EM.py.
The only change from the original EM is that the classifier
in the M-step for beta is an nn.Module optimised with gradient
descent instead of a logistic regression optimised with L-BFGS.

Everything else — E-step formula, M-step for T, convergence monitor,
the Dataset/EMResult dataclasses — is identical to T_estimation_EM.py.

Model
-----
    p(Y = k | X; theta)  =  softmax(MLP(X; theta))_k
    p(Y~= k | X; theta, T)  =  [p_Y(X; theta) @ T^T]_k

Contamination variants
----------------------
    "uniform"  :  T(eps)[k,l] = (1-eps)*I[k==l] + eps/K   (randomised response)
    "general"  :  T is a free column-stochastic matrix

EM iteration
------------
    E-step : gamma[i,l] = P(Y_i=l | Y~_i, X_i; theta^t, T^t)
    M-step : (a) update theta by minimising soft-CE loss against gamma
             (b) update T by closed-form update (identical to T_estimation_EM.py)
"""


# ---------------------------------------------------------------------------
# Contamination matrix helpers
# ---------------------------------------------------------------------------

def rr_matrix_torch(eps: torch.Tensor, K: int,
                    device: torch.device) -> torch.Tensor:
    eye  = torch.eye(K,  device=device)
    ones = torch.ones(K, K, device=device)
    return (1.0 - eps) * eye + (eps / K) * ones

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
# Classifier: any nn.Module that maps X -> logits over K classes
# ---------------------------------------------------------------------------

class MLP(nn.Module):
    """
    Default backbone: X -> logits (shape [batch, K]).
    Pass any other nn.Module to run_em() via the `classifier` argument.
    """
    def __init__(self, input_dim: int, K: int,
                 hidden_dims: list[int] = [128, 64]):
        super().__init__()
        layers = []
        prev = input_dim
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.ReLU()]
            prev = h
        layers.append(nn.Linear(prev, K))
        self.net = nn.Sequential(*layers)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        return self.net(X)
    
# ---------------------------------------------------------------------------
# Contamination parameter: thin wrapper so T has a unified interface
# ---------------------------------------------------------------------------

class RRContamination(nn.Module):
    """Randomised-response T, parametrised by a single scalar epsilon."""

    def __init__(self, K: int, eps_init: float = 0.3):
        super().__init__()
        self.K = K
        self.logit_eps = nn.Parameter(
            torch.tensor(eps_init).clamp(1e-6, 1 - 1e-6).logit()
        )

    @property
    def epsilon(self) -> torch.Tensor:
        return torch.sigmoid(self.logit_eps)

    def matrix(self) -> torch.Tensor:
        return rr_matrix_torch(self.epsilon, self.K, self.logit_eps.device)

    @torch.no_grad()
    def closed_form_update(self, gamma: torch.Tensor, Y_tilde: torch.Tensor):
        """
        Exact M-step for epsilon (eq. M-eps in the notes):
            eps* = K/(K-1) * (1 - mean_i gamma[i, Y~_i])
        """
        K = self.K
        diag = gamma[torch.arange(len(Y_tilde)), Y_tilde].mean().item()
        eps_new = float(np.clip((K / (K - 1)) * (1.0 - diag), 1e-6, 1 - 1e-6))
        self.logit_eps.data = torch.tensor(eps_new).logit()

    def eps_value(self) -> float:
        return self.epsilon.item()


class GeneralContamination(nn.Module):
    """
    General column-stochastic T, parametrised via softmax reparametrisation:
        T[:, l] = softmax(Psi[:, l])
    """

    def __init__(self, K: int):
        super().__init__()
        self.K = K
        self.Psi = nn.Parameter(torch.zeros(K, K))

    def matrix(self) -> torch.Tensor:
        return F.softmax(self.Psi, dim=0)   # normalise each column

    @torch.no_grad()
    def closed_form_update(self, gamma: torch.Tensor, Y_tilde: torch.Tensor,
                           floor: float = 1e-8):
        """
        Exact M-step for general T (remark in the notes):
            N[k, l]   = sum_{i: Y~_i = k} gamma[i, l]
            T[:, l]   = N[:, l] / sum_k N[k, l]
        We set Psi = log(T_new) so softmax(Psi) recovers T_new exactly.
        """
        K = self.K
        N = torch.zeros(K, K, device=gamma.device)
        N.index_add_(0, Y_tilde, gamma)
        col_sums = N.sum(dim=0, keepdim=True).clamp(min=1e-12)
        T_new = (N / col_sums).clamp(min=floor)
        T_new = T_new / T_new.sum(dim=0, keepdim=True)
        self.Psi.data = torch.log(T_new)

    def eps_value(self) -> float:
        """Approximate epsilon as average off-diagonal mass."""
        T = self.matrix().detach()
        K = self.K
        diag_mean = T.diag().mean().item()
        return float(np.clip((K / (K - 1)) * (1.0 - diag_mean), 0.0, 1.0))


# ---------------------------------------------------------------------------
# E-step
# ---------------------------------------------------------------------------

@torch.no_grad()
def e_step_nn(
    classifier: nn.Module,
    cont_layer:  nn.Module,        # RRContamination or GeneralContamination
    X_cont:      torch.Tensor,     # (n_c, d)
    Y_tilde:     torch.Tensor,     # (n_c,)  long
) -> torch.Tensor:
    """
    Compute responsibilities (same formula as T_estimation_EM.py):

        gamma[i, l] = T[Y~_i, l] * p_l(X_i; theta)
                      ---------------------------------
                      sum_{l'} T[Y~_i, l'] * p_{l'}(X_i; theta)

    Returns gamma of shape (n_c, K), detached (no gradient).
    """
    classifier.eval()
    log_p   = F.log_softmax(classifier(X_cont), dim=1)   # (n_c, K)
    T       = cont_layer.matrix()                          # (K, K)
    log_T   = torch.log(T[Y_tilde, :] + 1e-12)           # (n_c, K)

    log_gamma = log_p + log_T
    log_gamma = log_gamma - log_gamma.max(dim=1, keepdim=True).values
    gamma     = torch.exp(log_gamma)
    gamma     = gamma / gamma.sum(dim=1, keepdim=True)
    return gamma.detach()


# ---------------------------------------------------------------------------
# M-step for the classifier
# ---------------------------------------------------------------------------

def _classifier_loss(
    classifier: nn.Module,
    X_clean:    torch.Tensor,   # (n_cl, d)
    Y_clean:    torch.Tensor,   # (n_cl,)  hard labels
    X_cont:     torch.Tensor,   # (n_c, d)
    gamma:      torch.Tensor,   # (n_c, K)  soft targets, frozen
) -> torch.Tensor:
    """
    Q-function terms that depend on theta (eq. M-beta in the notes):

        - sum_{i in O}  log p_{Y_i}(X_i; theta)             [standard CE]
        - sum_{i in C}  sum_l gamma[i,l] log p_l(X_i; theta) [soft CE]

    The contamination layer does NOT appear here.
    gamma enters only as a fixed weighting — it is not differentiated through.
    """
    loss = torch.tensor(0.0, device=X_cont.device
                        if len(X_cont) else X_clean.device)

    if len(X_clean):
        loss = loss + F.cross_entropy(classifier(X_clean), Y_clean)

    if len(X_cont):
        log_p    = F.log_softmax(classifier(X_cont), dim=1)  # (n_c, K)
        soft_ce  = -(gamma * log_p).sum(dim=1).mean()
        loss     = loss + soft_ce

    return loss


def m_step_classifier(
    classifier:       nn.Module,
    X_clean:          torch.Tensor,
    Y_clean:          torch.Tensor,
    X_cont:           torch.Tensor,
    gamma:            torch.Tensor,
    optimizer:        torch.optim.Optimizer,
    n_steps:          int = 50,
    batch_size:       int = 256,
) -> float:

    #Maximise Q w.r.t. theta by running n_steps gradient steps.

    #Supports both Adam-style optimizers and L-BFGS.

    #L-BFGS note: because L-BFGS may evaluate the loss multiple times per
    #step (line search), it requires a closure. When L-BFGS is detected the
    #full dataset is used (no mini-batching) and n_steps maps to L-BFGS
    #iterations, exactly mirroring the original EM M-step for beta.

    #Returns the final loss value for monitoring.

    classifier.train()
    last_loss = float("nan")

    is_lbfgs = isinstance(optimizer, torch.optim.LBFGS)

    if is_lbfgs:
        # L-BFGS operates on the full dataset via a closure.
        # n_steps is handled internally by max_iter set on the optimizer.
        def closure():
            optimizer.zero_grad()
            loss = _classifier_loss(classifier, X_clean, Y_clean,
                                    X_cont, gamma)
            loss.backward()
            return loss

        loss_tensor = optimizer.step(closure)
        last_loss = loss_tensor.item()

    else:
        # Adam / SGD: mini-batch loop for n_steps steps.
        cont_loader = DataLoader(
            TensorDataset(X_cont, gamma), batch_size=batch_size, shuffle=True
        )
        step = 0
        while step < n_steps:
            for X_co_batch, gamma_batch in cont_loader:
                loss = _classifier_loss(classifier, X_clean, Y_clean,
                                        X_co_batch, gamma_batch)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                last_loss = loss.item()
                step += 1
                if step >= n_steps:
                    break

    return last_loss


# ---------------------------------------------------------------------------
# Observed log-likelihood monitor
# ---------------------------------------------------------------------------

@torch.no_grad()
def observed_log_likelihood_nn(
    classifier: nn.Module,
    cont_layer: nn.Module,
    X_clean:    torch.Tensor,
    Y_clean:    torch.Tensor,
    X_cont:     torch.Tensor,
    Y_tilde:    torch.Tensor,
) -> float:
    """
    Observed-data log-likelihood:

        sum_{i in O} log p_{Y_i}(X_i; theta)
      + sum_{i in C} log sum_l T[Y~_i, l] * p_l(X_i; theta)
    """
    classifier.eval()
    ll = 0.0

    if len(X_clean):
        log_p = F.log_softmax(classifier(X_clean), dim=1)
        ll   += log_p[torch.arange(len(Y_clean)), Y_clean].sum().item()

    if len(X_cont):
        p        = F.softmax(classifier(X_cont), dim=1)   # (n_c, K)
        T        = cont_layer.matrix()                      # (K, K)
        marginal = (T[Y_tilde, :] * p).sum(dim=1).clamp(min=1e-12)
        ll      += marginal.log().sum().item()

    return ll


# ---------------------------------------------------------------------------
# Result container  (same interface as T_estimation_EM.py)
# ---------------------------------------------------------------------------

@dataclass
class EMResult:
    classifier:      nn.Module
    cont_layer:      nn.Module
    log_likelihoods: list[float]
    n_iter:          int
    converged:       bool

    @property
    def T(self) -> np.ndarray:
        """Estimated contamination matrix as numpy array."""
        return self.cont_layer.matrix().detach().cpu().numpy()

    @property
    def eps(self) -> float:
        return self.cont_layer.eps_value()

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        Xt = torch.tensor(X, dtype=torch.float32)
        with torch.no_grad():
            return F.softmax(self.classifier(Xt), dim=1).cpu().numpy()

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.predict_proba(X).argmax(axis=1)


# ---------------------------------------------------------------------------
# Full EM loop
# ---------------------------------------------------------------------------

def run_em_nn(
    data:                Dataset,
    contamination_model: Literal["uniform", "general"] = "uniform",
    classifier:          Optional[nn.Module] = None,
    hidden_dims:         list[int] = [128, 64],
    eps_init:            float = 0,
    max_iter:            int = 100,
    tol:                 float = 1e-6,
    n_steps_per_iter:    int = 50,
    batch_size:          int = 256,
    lr:                  float = 1e-3,
    optimizer_name:      str = "adam",
    device:              str = "cpu",
    verbose:             bool = True,
) -> EMResult:
    """
    EM algorithm with an nn.Module classifier.

    Parameters
    ----------
    data                : Dataset (same as T_estimation_EM.py)
    contamination_model : "uniform" (RR) or "general"
    classifier          : any nn.Module mapping X -> logits of shape (batch, K).
                          If None, a default MLP with `hidden_dims` is used.
    hidden_dims         : MLP hidden layer sizes (ignored if classifier is given)
    eps_init            : initial contamination rate
    max_iter            : maximum EM iterations
    tol                 : convergence threshold on |Delta log-lik|
    n_steps_per_iter    : gradient steps on the classifier per M-step
                          (more steps → closer to exact EM M-step)
    batch_size          : mini-batch size for the classifier M-step
    lr                  : learning rate for the classifier optimizer
    device              : torch device string
    verbose             : print iteration log

    Returns
    -------
    """
    dev = torch.device(device)
    n, d = data.X.shape
    K    = data.K

    # --- Build or validate classifier ---
    if classifier is None:
        classifier = MLP(d, K, hidden_dims)
    classifier = classifier.to(dev)

    # --- Build contamination layer ---
    if contamination_model == "uniform":
        cont_layer = RRContamination(K, eps_init).to(dev)
    elif contamination_model == "general":
        cont_layer = GeneralContamination(K).to(dev)
    else:
        raise ValueError(f"Unknown contamination model: {contamination_model!r}")

    # --- Move data to torch tensors on device ---
    X_all   = torch.tensor(data.X,     dtype=torch.float32, device=dev)
    Y_all   = torch.tensor(data.Y_obs, dtype=torch.long,    device=dev)

    clean_mask = torch.tensor(data.clean_mask, device=dev)
    cont_mask  = torch.tensor(data.cont_mask,  device=dev)

    X_clean = X_all[clean_mask]
    Y_clean = Y_all[clean_mask]
    X_cont  = X_all[cont_mask]
    Y_tilde = Y_all[cont_mask]

    # Classifier optimizer (contamination layer updated separately / analytically)
    if optimizer_name == "lbfgs":
        optimizer = torch.optim.LBFGS(
            classifier.parameters(), lr=lr,
            max_iter=n_steps_per_iter, line_search_fn="strong_wolfe"
        )
    else:
        optimizer = torch.optim.Adam(classifier.parameters(), lr=lr)

    log_likelihoods: list[float] = []
    converged = False

    if verbose:
        header = f"{'Iter':>5}  {'Log-lik':>14}  {'Delta':>12}  {'eps':>8}"
        print(header)
        print("-" * (len(header) + 2))

    for t in range(max_iter):

        # ------------------------------------------------------------------
        # E-step: compute responsibilities for contaminated samples
        # ------------------------------------------------------------------
        gamma = e_step_nn(classifier, cont_layer, X_cont, Y_tilde)  # (n_c, K)

        # ------------------------------------------------------------------
        # M-step (a): update classifier with contamination layer frozen
        # ------------------------------------------------------------------
        m_step_classifier(
            classifier, X_clean, Y_clean, X_cont, gamma,
            optimizer, n_steps=n_steps_per_iter, batch_size=batch_size,
        )

        # ------------------------------------------------------------------
        # M-step (b): update T with closed-form update (exact EM M-step)
        # ------------------------------------------------------------------
        cont_layer.closed_form_update(gamma, Y_tilde)

        # ------------------------------------------------------------------
        # Monitor observed log-likelihood
        # ------------------------------------------------------------------
        ll = observed_log_likelihood_nn(
            classifier, cont_layer, X_clean, Y_clean, X_cont, Y_tilde
        )
        log_likelihoods.append(ll)

        if verbose:
            delta = ll - log_likelihoods[-2] if t > 0 else float("nan")
            print(f"{t+1:>5}  {ll:>14.4f}  {delta:>12.6f}  "
                  f"{cont_layer.eps_value():>8.4f}")

        # ------------------------------------------------------------------
        # Convergence check
        # ------------------------------------------------------------------
        if t > 0 and abs(ll - log_likelihoods[-2]) < tol:
            if verbose:
                print(f"\nConverged at iteration {t + 1}.")
            converged = True
            break
    else:
        if verbose:
            warnings.warn(f"EM did not converge in {max_iter} iterations.")

    return EMResult(
        classifier=classifier,
        cont_layer=cont_layer,
        log_likelihoods=log_likelihoods,
        n_iter=len(log_likelihoods),
        converged=converged,
    )


# ---------------------------------------------------------------------------
# Prediction helpers
# ---------------------------------------------------------------------------

def predict_proba(classifier: nn.Module, X: np.ndarray) -> np.ndarray:
    Xt = torch.tensor(X, dtype=torch.float32)
    with torch.no_grad():
        return F.softmax(classifier(Xt), dim=1).cpu().numpy()

def predict(classifier: nn.Module, X: np.ndarray) -> np.ndarray:
    return predict_proba(classifier, X).argmax(axis=1)
