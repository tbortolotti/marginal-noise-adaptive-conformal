from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset

def compute_ll_cont(model, X_feat_torch, Yt_torch):
    """Marginal log-likelihood on contaminated observations."""
    model.eval()
    with torch.no_grad():
        logits = model.backbone(X_feat_torch)
        p = F.softmax(logits, dim=-1)
        T = model.contamination.contamination_matrix()
        marginal = (T[Yt_torch, :] * p).sum(dim=1).clamp(min=1e-12)
        return marginal.log().mean().item()


class RandomizedResponseLayer(nn.Module):
    """
    Differentiable layer implementing the Randomized Response contamination.

    Given a distribution p_Y (shape: [batch, K]), computes the contaminated
    distribution p_Ytilde = T(epsilon) @ p_Y, where:

        T[j, k] = epsilon * I[j==k] + (1 - epsilon) / K

    epsilon is a learnable scalar, constrained to [0, 1] via sigmoid.
    """

    def __init__(self, K: int, epsilon_init: float = 0.5):
        super().__init__()
        self.K = K
        # Parametrize epsilon via logit so that sigma(logit_eps) in (0,1)
        logit_init = torch.tensor(epsilon_init).clamp(1e-6, 1 - 1e-6).logit()
        self.logit_epsilon = nn.Parameter(logit_init)

    @property
    def epsilon(self) -> torch.Tensor:
        """Returns epsilon in (0, 1)."""
        return torch.sigmoid(self.logit_epsilon)

    def contamination_matrix(self) -> torch.Tensor:
        """
        Builds T(epsilon) of shape [K, K].
        T[j, k] = (1-epsilon) * I[j==k] + epsilon / K
        """
        eps = self.epsilon
        ones = torch.ones(self.K, self.K, device=self.logit_epsilon.device)
        eye = torch.eye(self.K, device=self.logit_epsilon.device)
        T = eps / self.K * ones + (1.0 - eps) * eye

        #eps = self.epsilon
        #T = (1.0 - eps) / self.K * torch.ones(self.K, self.K,
        #                                       device=self.logit_epsilon.device)
        #T = T + eps * torch.eye(self.K, device=self.logit_epsilon.device)
        #return T
        
        return T

    def forward(self, p_Y: torch.Tensor) -> torch.Tensor:
        """
        Args:
            p_Y: clean label distribution, shape [batch, K]
        Returns:
            p_Ytilde: contaminated distribution, shape [batch, K]
        """
        T = self.contamination_matrix()          # [K, K]
        p_Ytilde = p_Y @ T.t()                  # [batch, K] @ [K, K] = [batch, K]
        return p_Ytilde
    

    # ------------------------------------------------------------------
    # Closed-form M-step for the RR model (mirrors EM eq. M-eps)
    # ------------------------------------------------------------------
    @torch.no_grad()
    def closed_form_update(self, gamma: torch.Tensor,
                           Y_tilde: torch.Tensor,
                           eps_bounds: tuple = (1e-6, 1.0 - 1e-6)):
        """
        gamma   : (n_c, K)  responsibilities computed in the E-step
        Y_tilde : (n_c,)    contaminated labels (long, 0-indexed)

        Closed-form update (from EM M-step for epsilon):
            eps* = K/(K-1) * (1 - mean_i gamma[i, Y~_i])
        """
        K = self.K
        diag_weight = gamma[torch.arange(len(Y_tilde)), Y_tilde].mean().item()
        eps_new = (K / (K - 1)) * (1.0 - diag_weight)
        eps_new = float(np.clip(eps_new, *eps_bounds))
        # Back-solve: eps_new = sigmoid(phi)  =>  phi = logit(eps_new)
        self.logit_epsilon.data = torch.tensor(eps_new).logit()
    

class GeneralContaminationLayer(nn.Module):

    def __init__(self, K: int,
                 epsilon_init: float = 0.0,
                 floor: float = 1e-6):
        super().__init__()
        self.K = K

        # Build RR-style initialization matrix
        eps = float(np.clip(epsilon_init, floor, 1.0 - floor))
        T_init = eps / K * torch.ones(K, K) + (1.0 - eps) * torch.eye(K)

        # Free parameters: one K-vector per column of T
        # Shape [K, K]: entry [k, l] = psi_{l,k}
        Psi_init = torch.log(T_init)
        self.Psi = nn.Parameter(Psi_init)

    def contamination_matrix(self) -> torch.Tensor:
        # Apply softmax column by column
        # softmax over dim=0 normalizes each column
        return F.softmax(self.Psi, dim=0)   # shape [K, K]

    def forward(self, p_Y: torch.Tensor) -> torch.Tensor:
        T = self.contamination_matrix()
        return p_Y @ T.t()                  # [batch, K]
    

    # ------------------------------------------------------------------
    # Closed-form M-step for general T (mirrors EM remark on general T)
    # ------------------------------------------------------------------
    @torch.no_grad()
    def closed_form_update(self, gamma: torch.Tensor,
                           Y_tilde: torch.Tensor,
                           floor: float = 1e-8):
        """
        Aggregated responsibility matrix:
            N[k, l] = sum_{i : Y~_i = k}  gamma[i, l]

        Closed-form column update:
            T[:, l] = N[:, l] / sum_k N[k, l]

        We set Psi so that softmax(Psi, dim=0) recovers T_new,
        i.e. Psi[:, l] = log(T_new[:, l])  (up to a constant).
        """
        K = self.K
        N = torch.zeros(K, K, device=gamma.device)
        N.index_add_(0, Y_tilde, gamma)       # N[Y~_i, :] += gamma[i, :]

        col_sums = N.sum(dim=0, keepdim=True).clamp(min=1e-12)
        T_new = (N / col_sums).clamp(min=floor)
        T_new = T_new / T_new.sum(dim=0, keepdim=True)  # renormalise

        # Recover Psi such that softmax(Psi, dim=0) = T_new
        self.Psi.data = torch.log(T_new)


# ---------------------------------------------------------------------------
# Backbone classifier
# ---------------------------------------------------------------------------
class ClassifierBackbone(nn.Module):
    """
    MLP backbone mapping X -> logits over Y.
    Architecture: X -> Linear -> ReLU -> ... -> Linear -> logits (K)
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
        """Returns logits of shape [batch, K]."""
        return self.net(X)

# ---------------------------------------------------------------------------
# Full model
# ---------------------------------------------------------------------------
class NoisyLabelNet(nn.Module):
    """
    Full model combining backbone + contamination layer.

    Forward pass:
        X -> (logits_Y, logits_Ytilde)
    """

    def __init__(self, input_dim: int = None,
                 K: int = 2,
                 hidden_dims: list[int] = [128, 64],
                 contamination_model_: str = "uniform",
                 epsilon_init: float = 0):
        super().__init__()
        self.K = K

        self.backbone = ClassifierBackbone(input_dim, K, hidden_dims)

        if contamination_model_=="uniform":
            self.contamination = RandomizedResponseLayer(K, epsilon_init)
        elif contamination_model_=="general":
            self.contamination = GeneralContaminationLayer(K, epsilon_init)

    def forward(self, X: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            X     : features, shape [batch, input_dim]
        Returns:
            logits_Y      : logits for Y,  shape [batch, K]
            logits_Ytilde : logits for Ỹ, shape [batch, K]
        """
        # --- Backbone: X -> logits_Y -> p_Y ---
        logits_Y = self.backbone(X)                       # [batch, K]
        p_Y = F.softmax(logits_Y, dim=-1)                 # [batch, K]
        p_Ytilde = self.contamination(p_Y)

        # Convert back to log-space for loss computation
        logits_Ytilde = torch.log(p_Ytilde + 1e-8)        # [batch, K]

        return logits_Y, logits_Ytilde

    @property
    def epsilon(self) -> float:
        return self.contamination.epsilon.item()


# ---------------------------------------------------------------------------
# Loss function
# ---------------------------------------------------------------------------
def noisy_label_loss(logits_Y: torch.Tensor,
                     logits_Ytilde: torch.Tensor,
                     obs_labels: torch.Tensor,
                     I: torch.Tensor,
                     loss_type: str = "equal") -> torch.Tensor:
    """
    Selective cross-entropy loss:
        - Use logits_Y     for samples where I == 1 (clean label observed)
        - Use logits_Ytilde for samples where I == 0 (contaminated label observed)

    Args:
        logits_Y      : shape [batch, K]
        logits_Ytilde : shape [batch, K]  (already in log-space)
        obs_labels    : shape [batch], integer class indices
        I             : shape [batch], binary indicator (1 = clean, 0 = noisy)
    Returns:
        scalar loss
    """

    if loss_type=="equal":
        I = I.bool()
        loss = torch.zeros(logits_Y.shape[0], device=logits_Y.device)

        # Cross-entropy for clean observations
        if I.any():
            loss[I] = F.cross_entropy(logits_Y[I], obs_labels[I], reduction='none')

        # Cross-entropy for contaminated observations
        # logits_Ytilde is already log-softmax-like, use nll_loss
        if (~I).any():
            loss[~I] = F.nll_loss(logits_Ytilde[~I], obs_labels[~I], reduction='none')

        loss_ = loss.mean()
    elif loss_type=="weighted":
        I = I.bool()
        losses = []

        # Cross-entropy for clean observations
        if I.any():
            loss_clean = F.cross_entropy(logits_Y[I], obs_labels[I], reduction='mean')
            losses.append(loss_clean)

        # Cross-entropy for contaminated observations
        if (~I).any():
            loss_noisy = F.nll_loss(logits_Ytilde[~I], obs_labels[~I], reduction='none')
            losses.append(loss_noisy.mean())
        
        loss_ = sum(losses) / len(losses)
    elif loss_type=="upweighted":
        I = I.bool()
        pi_clean_ = torch.round(torch.sum(I==1)/len(I), decimals=2)
        losses = []

        # Cross-entropy for clean observations
        if I.any():
            loss_clean = F.cross_entropy(logits_Y[I], obs_labels[I], reduction='mean')
            loss_clean_ = 1/pi_clean_*loss_clean
            losses.append(loss_clean_)

        # Cross-entropy for contaminated observations
        if (~I).any():
            loss_noisy = F.nll_loss(logits_Ytilde[~I], obs_labels[~I], reduction='none')
            loss_noisy_ = 1/(1-pi_clean_) * loss_noisy.mean()
            losses.append(loss_noisy_)
        
        loss_ = sum(losses) / len(losses)

    return loss_

def contamination_regularization(model, lambda_reg=0.1):
    T_current = model.contamination.contamination_matrix()

    # Maximize log|det(T)| — directly encourages invertibility
    # Use SVD for numerical stability
    sign, logabsdet = torch.linalg.slogdet(T_current)
    # We want det > 0 and large, so penalize -log|det|
    reg = -logabsdet  # minimizing this maximizes |det(T)|
        
    return lambda_reg * reg

# ---------------------------------------------------------------------------
# Training loop
# ---------------------------------------------------------------------------
def train(model: NoisyLabelNet,
          X: torch.Tensor,
          obs_labels: torch.Tensor,
          I: torch.Tensor,
          n_epochs: int = 100,
          batch_size: int = 256,
          lr: float = 1e-3,
          device: str = "cpu",
          loss_type: str = "equal",
          verbose: bool = False) -> list[dict]:

    model = model.to(device)
    X, obs_labels, I = X.to(device), obs_labels.to(device), I.to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)
    dataset = TensorDataset(X, obs_labels, I)
    loader  = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    history = []
    for epoch in range(n_epochs):
        model.train()
        total_loss = 0.0

        for X_batch, labels_batch, I_batch in loader:
            logits_Y, logits_Ytilde = model(X_batch)
            loss = noisy_label_loss(logits_Y, logits_Ytilde, labels_batch, I_batch, loss_type)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item() * X_batch.shape[0]

        avg_loss = total_loss / len(dataset)
        #eps_val  = model.epsilon

        if verbose:
            if (epoch + 1) % 10 == 0:
                #print(f"Epoch {epoch+1:4d}/{n_epochs} | loss={avg_loss:.4f} | ε={eps_val:.4f}")
                print(f"Epoch {epoch+1:4d}/{n_epochs} | loss={avg_loss:.4f}")

        #history.append({"epoch": epoch + 1, "loss": avg_loss, "epsilon": eps_val})
        history.append({"epoch": epoch + 1, "loss": avg_loss})

    return history


# ---------------------------------------------------------------------------
# Alternate training loop
# 1. Freeze contamination, run several gradient steps on backbone
# 2. Freeze backbone, run several gradient steps on contamination
# ---------------------------------------------------------------------------

def train_alternate(model: NoisyLabelNet,
                    X: torch.Tensor,
                    obs_labels: torch.Tensor,
                    I: torch.Tensor,
                    n_epochs: int = 100,
                    n_grad_steps: int = 50,
                    batch_size: int = 256,
                    lr: float = 1e-3,
                    lambda_reg: float = 0,
                    device: str = "cpu",
                    loss_type: str = "equal",
                    verbose: bool = False) -> list[dict]:

    model = model.to(device)
    X, obs_labels, I = X.to(device), obs_labels.to(device), I.to(device)

    cont_mask  = (I == 0)
    X_cont  = X[cont_mask]
    I_cont = I[cont_mask]
    obs_labels_cont = obs_labels[cont_mask]

    optimizer_backbone      = torch.optim.AdamW(model.backbone.parameters(), lr=lr)
    optimizer_contamination = torch.optim.AdamW(model.contamination.parameters(), lr=lr)

    dataset = TensorDataset(X, obs_labels, I)
    history = []

    for epoch in range(n_epochs):
        model.train()
        total_loss_1 = 0.0
        total_loss_2 = 0.0

        # ------------------------------------------------------------------
        # Phase 1: freeze contamination, update backbone
        # ------------------------------------------------------------------
        for p in model.backbone.parameters():
            p.requires_grad_(True)
        for p in model.contamination.parameters():
            p.requires_grad_(False)

        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
        step = 0
        while step < n_grad_steps:
            for X_batch, labels_batch, I_batch in loader:
                logits_Y, logits_Ytilde = model(X_batch)
                loss = noisy_label_loss(logits_Y, logits_Ytilde,
                                        labels_batch, I_batch, loss_type)
                optimizer_backbone.zero_grad()
                loss.backward()
                optimizer_backbone.step()
                total_loss_1 += loss.item()
                step += 1
                if step >= n_grad_steps:
                    break

        # ------------------------------------------------------------------
        # Phase 2: freeze backbone, update contamination
        # ------------------------------------------------------------------
        for p in model.backbone.parameters():
            p.requires_grad_(False)
        for p in model.contamination.parameters():
            p.requires_grad_(True)

        cont_dataset = TensorDataset(X_cont, obs_labels_cont, I_cont)
        cont_loader  = DataLoader(cont_dataset, batch_size=batch_size, shuffle=True)

        step = 0
        for _ in range(n_grad_steps):
            for X_batch, labels_batch, I_batch in cont_loader:
                logits_Y, logits_Ytilde = model(X_batch)
                loss = noisy_label_loss(logits_Y, logits_Ytilde,
                                        labels_batch, I_batch, loss_type)
                loss = loss + contamination_regularization(model, lambda_reg=lambda_reg)
                optimizer_contamination.zero_grad()
                loss.backward()
                optimizer_contamination.step()
                step += 1
                if step >= n_grad_steps:
                    break

        # Restore all gradients
        for p in model.parameters():
            p.requires_grad_(True)

        if verbose and (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1:4d}/{n_epochs} | "
                  f"loss_backbone={total_loss_1:.4f} | "
                  f"loss_cont={total_loss_2:.4f}")

        history.append({"epoch": epoch + 1,
                        "loss_backbone": total_loss_1,
                        "loss_cont": total_loss_2})

    return history


