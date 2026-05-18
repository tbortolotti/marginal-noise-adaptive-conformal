import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
#import torchvision.models as models
import sys
sys.path.append("/home1/tb_214/code/PyTorch_CIFAR10")
from cifar10_models.resnet import resnet18 as cifar_resnet18

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

"""
class ResNetBackbone(nn.Module):
    def __init__(self, K: int, freeze_features: bool = False):
        super().__init__()
        resnet = models.resnet18(weights=models.ResNet18_Weights.DEFAULT)
        
        # Replace the final FC layer to output K logits
        in_features = resnet.fc.in_features  # 512 for ResNet-18
        resnet.fc = nn.Linear(in_features, K)
        
        self.net = resnet
        
        if freeze_features:
            # Freeze everything except the final FC layer
            for name, param in self.net.named_parameters():
                if name not in ("fc.weight", "fc.bias"):
                    param.requires_grad_(False)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        return self.net(X)
"""

class ResNetBackbone(nn.Module):
    def __init__(self, K: int, freeze_features: bool = False):
        super().__init__()
        self.net = cifar_resnet18(pretrained=True)
        
        if freeze_features:
            for name, param in self.net.named_parameters():
                if "linear" not in name:
                    param.requires_grad_(False)


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
                 backbone_model_: str = "MLP",
                 hidden_dims: list[int] = [128, 64],
                 freeze_features: bool = False,
                 contamination_model_: str = "uniform",
                 epsilon_init: float = 0):
        super().__init__()
        self.K = K

        if backbone_model_=="resnet":
            self.backbone = ResNetBackbone(K, freeze_features=freeze_features)
        else:
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



# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# EM-style Alternating Training of the NN
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# E-step  (computes soft responsibilities — no gradient)
# ---------------------------------------------------------------------------

@torch.no_grad()
def e_step_nn(
    model: NoisyLabelNet,
    X_cont: torch.Tensor,    # (n_c, d)
    Y_tilde: torch.Tensor,   # (n_c,)  contaminated labels, long
) -> torch.Tensor:
    """
    Compute responsibilities in PyTorch, mirroring the EM E-step:

        gamma[i, l]  proportional to  T[Y~_i, l] * p_l(X_i; theta)

    Returns gamma of shape (n_c, K), no gradient.
    """
    model.eval()
    logits_Y = model.backbone(X_cont)          # (n_c, K)
    log_p    = F.log_softmax(logits_Y, dim=1)  # (n_c, K)

    T        = model.contamination.contamination_matrix()   # (K, K)
    log_T    = torch.log(T[Y_tilde, :] + 1e-12)            # (n_c, K)

    log_gamma = log_p + log_T
    log_gamma = log_gamma - log_gamma.max(dim=1, keepdim=True).values
    gamma     = torch.exp(log_gamma)
    gamma     = gamma / gamma.sum(dim=1, keepdim=True)
    return gamma   # (n_c, K), detached


# ---------------------------------------------------------------------------
# M-step for backbone: soft CE loss with frozen contamination
# ---------------------------------------------------------------------------

def backbone_loss(
    model: NoisyLabelNet,
    X_batch: torch.Tensor,
    labels_batch: torch.Tensor,
    I_batch: torch.Tensor,
    gamma_batch: torch.Tensor,   # same size as batch, zeros for clean
) -> torch.Tensor:

    logits = model.backbone(X_batch)
    log_p  = F.log_softmax(logits, dim=1)

    I_bool = I_batch.bool()

    loss = torch.tensor(0.0, device=X_batch.device)

    # Clean part (standard CE)
    if I_bool.any():
        loss_clean = F.cross_entropy(
            logits[I_bool],
            labels_batch[I_bool],
            reduction='mean'
        )
        loss = loss + loss_clean

    # Contaminated part (soft CE)
    if (~I_bool).any():
        soft_ce = -(gamma_batch[~I_bool] * log_p[~I_bool]).sum(dim=1).mean()
        loss = loss + soft_ce

    return loss


# ---------------------------------------------------------------------------
# EM-style alternating training loop
# ---------------------------------------------------------------------------

def train_em_style(
    model: NoisyLabelNet,
    X: torch.Tensor,
    obs_labels: torch.Tensor,
    I: torch.Tensor,
    n_epochs: int = 100,
    optimizer_name: str = "adam",
    # Number of gradient steps on the backbone per EM iteration.
    # Mimics the inner L-BFGS convergence of the EM M-step.
    n_backbone_steps: int = 50,
    batch_size: int = 256,
    lr_backbone: float = 1e-3,
    device: str = "cpu",
    verbose: bool = True,
) -> list[dict]:
    """
    EM-style alternating training for NoisyLabelNet.

    Each "epoch" corresponds to one EM iteration:

        1. E-step  : compute soft responsibilities gamma for contaminated
                     samples using current (theta, T). No gradient.
        2. M-step backbone : freeze contamination, take n_backbone_steps
                     gradient steps minimising backbone_loss(theta; gamma).
                     This optimises the Q-function w.r.t. theta exactly as
                     the EM M-step for beta does.
        3. M-step T: freeze backbone.
                     - If use_closed_form_cont=True: apply the closed-form
                       update (identical to the EM M-step for epsilon/T).
                     - Otherwise: take n_grad_steps gradient steps on the
                       contamination layer parameters.

    Key difference from joint SGD
    --------------------------------
    In joint training, theta and phi are updated simultaneously on every
    mini-batch, so the contamination layer never sees a coherent set of
    responsibilities. Here, gamma is recomputed once per epoch (E-step)
    and then held fixed while each parameter group is updated in turn,
    exactly mirroring the EM's alternating maximisation.
    """
    model = model.to(device)
    X, obs_labels, I = X.to(device), obs_labels.to(device), I.to(device)

    clean_mask = (I == 1)
    cont_mask  = (I == 0)

    X_clean = X[clean_mask]
    Y_clean = obs_labels[clean_mask]
    X_cont  = X[cont_mask]
    Y_tilde = obs_labels[cont_mask]

    # Separate optimizers — only parameters of the respective module are updated
    if optimizer_name == "adam":
        optimizer_backbone = torch.optim.AdamW(
            model.backbone.parameters(), lr=lr_backbone
        )
    elif optimizer_name == "lbfgs":
        optimizer_backbone = torch.optim.LBFGS(
            model.backbone.parameters(),
            lr=lr_backbone,
            max_iter=n_backbone_steps,
            history_size=10,
            line_search_fn="strong_wolfe",
        )
    else:
        raise ValueError("optimizer_name must be 'adam' or 'lbfgs'")

    history = []

    for epoch in range(n_epochs):

        # ==================================================================
        # E-step: compute responsibilities (no gradient, full dataset)
        # ==================================================================
        gamma = e_step_nn(model, X_cont, Y_tilde)   # (n_c, K), detached

        # Build full gamma tensor aligned with X
        gamma_full = torch.zeros(X.shape[0], model.K, device=X.device)
        gamma_full[cont_mask] = gamma

        # ==================================================================
        # M-step 1: update backbone with contamination frozen
        # ==================================================================
        for p in model.contamination.parameters():
            p.requires_grad_(False)
        for p in model.backbone.parameters():
            p.requires_grad_(True)

        model.train()
        #last_loss = float("nan")

        if optimizer_name == "adam":
            # Mini-batch
            dataset = TensorDataset(X, obs_labels, I, gamma_full)
            loader  = DataLoader(dataset, batch_size=batch_size, shuffle=True)

            step = 0
            for _ in range(n_backbone_steps):
                for X_b, y_b, I_b, g_b in loader:

                    loss = backbone_loss(model, X_b, y_b, I_b, g_b)

                    optimizer_backbone.zero_grad()
                    loss.backward()
                    optimizer_backbone.step()

                    step += 1
                    if step >= n_backbone_steps:
                        break
                if step >= n_backbone_steps:
                    break

        elif optimizer_name == "lbfgs":
            # Full batch
            def closure():
                optimizer_backbone.zero_grad()
                loss = backbone_loss(model, X, obs_labels, I, gamma_full)
                loss.backward()
                return loss

            for _ in range(n_backbone_steps):
                optimizer_backbone.step(closure)

        # ==================================================================
        # M-step 2: update contamination with backbone frozen
        # ==================================================================
        for p in model.backbone.parameters():
            p.requires_grad_(False)
        for p in model.contamination.parameters():
            p.requires_grad_(True)

        # Exact EM update — no gradient needed
        model.contamination.closed_form_update(gamma, Y_tilde)

        # Restore all grads for monitoring
        for p in model.parameters():
            p.requires_grad_(True)

        # ==================================================================
        # Diagnostics
        # ==================================================================
        with torch.no_grad():
            model.eval()
            # Observed log-likelihood (contaminated part)
            ll_clean = 0.0
            if len(X_clean):
                logits_cl = model.backbone(X_clean)
                ll_clean  = -F.cross_entropy(logits_cl, Y_clean).item()

            ll_cont = 0.0
            if len(X_cont):
                logits_co = model.backbone(X_cont)
                p_co      = F.softmax(logits_co, dim=-1)
                T         = model.contamination.contamination_matrix()
                marginal  = (T[Y_tilde, :] * p_co).sum(dim=1).clamp(min=1e-12)
                ll_cont   = marginal.log().mean().item()

            #eps_val = model.epsilon

        record = {
            "epoch":    epoch + 1,
            "ll_clean": ll_clean,
            "ll_cont":  ll_cont,
        }
        history.append(record)

        if verbose and (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1:4d}/{n_epochs} | "
                  f"ll_clean={ll_clean:.4f} | ll_cont={ll_cont:.4f}")

    return history


# ---------------------------------------------------------------------------
# NOTES (code that'll probably be discarded)
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
"""