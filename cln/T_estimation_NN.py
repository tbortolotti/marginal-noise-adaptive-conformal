import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset

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
        logit_init = torch.tensor(epsilon_init).logit()
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
    

class GeneralContaminationLayer(nn.Module):

    def __init__(self, K: int):
        super().__init__()
        self.K = K
        # Free parameters: one K-vector per column of T
        # Shape [K, K]: entry [k, l] = psi_{l,k}
        self.Psi = nn.Parameter(torch.zeros(K, K))

    def contamination_matrix(self) -> torch.Tensor:
        # Apply softmax column by column
        # softmax over dim=0 normalizes each column
        return F.softmax(self.Psi, dim=0)   # shape [K, K]

    def forward(self, p_Y: torch.Tensor) -> torch.Tensor:
        T = self.contamination_matrix()
        return p_Y @ T.t()                  # [batch, K]


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

    def __init__(self, input_dim: int, K: int,
                 hidden_dims: list[int] = [128, 64],
                 contamination_model_: str = "uniform",
                 epsilon_init: float = 0.5):
        super().__init__()
        self.K = K

        self.backbone = ClassifierBackbone(input_dim, K, hidden_dims)
        if contamination_model_=="uniform":
            self.contamination = RandomizedResponseLayer(K, epsilon_init)
        elif contamination_model_=="general":
            self.contamination = GeneralContaminationLayer(K)

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


# ---------------------------------------------------------------------------
# Training loop
# ---------------------------------------------------------------------------
"""
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
"""

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

    # Split indices into clean (I==1) and noisy (I==0)
    clean_idx = (I == 1).nonzero(as_tuple=True)[0]
    noisy_idx = (I == 0).nonzero(as_tuple=True)[0]

    # Separate optimizers for backbone and contamination layer
    optimizer_backbone     = torch.optim.AdamW(model.backbone.parameters(), lr=lr)
    optimizer_contamination = torch.optim.AdamW(model.contamination.parameters(), lr=lr)

    # How many steps per epoch: driven by the larger of the two splits
    steps_per_epoch = max(len(clean_idx), len(noisy_idx)) // batch_size + 1

    history = []
    for epoch in range(n_epochs):
        model.train()
        total_loss_clean = 0.0
        total_loss_noisy = 0.0

        # Shuffle split indices at the start of each epoch
        clean_perm = clean_idx[torch.randperm(len(clean_idx), device=device)]
        noisy_perm = noisy_idx[torch.randperm(len(noisy_idx), device=device)]

        for step in range(steps_per_epoch):
            # ----------------------------------------------------------------
            # Phase 1: clean batch → update backbone only
            # ----------------------------------------------------------------
            # Cycle over clean indices if exhausted
            start = (step * batch_size) % len(clean_perm)
            end   = start + batch_size
            if end > len(clean_perm):               # wrap around
                batch_idx = torch.cat([clean_perm[start:], clean_perm[:end - len(clean_perm)]])
            else:
                batch_idx = clean_perm[start:end]

            X_clean     = X[batch_idx]
            labels_clean = obs_labels[batch_idx]
            I_clean     = I[batch_idx]              # all ones, but kept for generality

            # Freeze contamination, unfreeze backbone
            for p in model.contamination.parameters():
                p.requires_grad_(False)
            for p in model.backbone.parameters():
                p.requires_grad_(True)

            logits_Y, logits_Ytilde = model(X_clean)
            loss_clean = noisy_label_loss(logits_Y, logits_Ytilde, labels_clean, I_clean, loss_type)

            optimizer_backbone.zero_grad()
            loss_clean.backward()
            optimizer_backbone.step()

            # ----------------------------------------------------------------
            # Phase 2: noisy batch → update contamination layer only
            # ----------------------------------------------------------------
            start = (step * batch_size) % len(noisy_perm)
            end   = start + batch_size
            if end > len(noisy_perm):
                batch_idx = torch.cat([noisy_perm[start:], noisy_perm[:end - len(noisy_perm)]])
            else:
                batch_idx = noisy_perm[start:end]

            X_noisy      = X[batch_idx]
            labels_noisy = obs_labels[batch_idx]
            I_noisy      = I[batch_idx]             # all zeros, but kept for generality

            # Freeze backbone, unfreeze contamination
            for p in model.backbone.parameters():
                p.requires_grad_(False)
            for p in model.contamination.parameters():
                p.requires_grad_(True)

            logits_Y, logits_Ytilde = model(X_noisy)
            loss_noisy = noisy_label_loss(logits_Y, logits_Ytilde, labels_noisy, I_noisy, loss_type)

            optimizer_contamination.zero_grad()
            loss_noisy.backward()
            optimizer_contamination.step()

            total_loss_clean += loss_clean.item()
            total_loss_noisy += loss_noisy.item()

        avg_loss_clean = total_loss_clean / steps_per_epoch
        avg_loss_noisy = total_loss_noisy / steps_per_epoch

        if verbose and (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1:4d}/{n_epochs} | "
                  f"loss_clean={avg_loss_clean:.4f} | "
                  f"loss_noisy={avg_loss_noisy:.4f}")

        history.append({"epoch": epoch + 1,
                        "loss_clean": avg_loss_clean,
                        "loss_noisy": avg_loss_noisy,
        })

    # Restore all gradients at the end
    for p in model.parameters():
        p.requires_grad_(True)

    return history