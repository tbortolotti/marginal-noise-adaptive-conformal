import numpy as np
from scipy.stats import norm
from sklearn.model_selection import train_test_split
from sklearn.datasets import make_classification
import pdb
from cln import contamination

def sigmoid(x):
    return(1/(1 + np.exp(-x)))

def sample_truncated_gaussian(n, R=1.0, A=None, mu=None, batch_size=50000, random_state=None):
    rng = np.random.default_rng(random_state)

    A = np.eye(2) if A is None else np.asarray(A, dtype=float)

    out = np.empty((n, 2), dtype=float)
    filled = 0

    while filled < n:
        m = min(batch_size, (n - filled) * 10)
        u = rng.normal(size=(m, 2))
        accept = (u[:, 0] * u[:, 0] + u[:, 1] * u[:, 1]) <= (R * R)

        accepted_u = u[accept]
        k = accepted_u.shape[0]
        if k == 0:
            continue

        take = min(k, n - filled)
        U_take = accepted_u[:take]

        X_take = (U_take @ A.T) + mu[filled:filled + take]
        out[filled:filled + take] = X_take
        filled += take

    return out

class DataModel:
    def __init__(self, K, p, random_state=None):
        self.K = K
        self.p = p
        self.rng = np.random.default_rng(seed=random_state)
        self.random_state = random_state
        
    def sample(self, n):
        X = self.sample_X(n)
        Y = self.sample_Y(X)
        return X, Y

    def set_seed(self, random_state):
        self.random_state = random_state
        self.rng = np.random.default_rng(seed=random_state)

    def estimate_rho(self, n_mc=10000):
        _, Y = self.sample(n_mc)
        rho = np.ones((self.K,))
        for k in range(self.K):
            rho[k] = np.mean(Y==k)
        return rho


class DataModel_1(DataModel):
    def __init__(self, K, p, signal=1, random_state=None):
        super().__init__(K, p, random_state=random_state)
        self.n_informative = np.maximum(1, int(self.p*0.5))
        self.class_sep = signal

    def sample(self, n):        
        X, Y = make_classification(n_samples=n, n_classes=self.K,
                                   n_features=self.p, n_informative=self.n_informative, 
                                   class_sep=self.class_sep, random_state=self.random_state)
        Y = Y.astype(np.int32)
        return X, Y
    
class DataModel_1_easy(DataModel):
    def __init__(self, K, p, signal=1, flipy=0.01, random_state=None):
        super().__init__(K, p, random_state=random_state)
        self.n_informative = 2
        self.n_redundant = p - self.n_informative
        self.n_clusters_per_class = 1
        self.class_sep = signal
        self.flipy = flipy

    def sample(self, n):        
        X, Y = make_classification(n_samples=n,
                                   n_classes=self.K,
                                   n_features=self.p, n_informative=self.n_informative,
                                   n_redundant=self.n_redundant,
                                   n_clusters_per_class=self.n_clusters_per_class,
                                   flip_y=self.flipy,
                                   class_sep=self.class_sep, random_state=self.random_state)
        Y = Y.astype(np.int32)
        return X, Y
    

class DataModel_2(DataModel):
    def __init__(self, K, p, signal=1, random_state=None):
        super().__init__(K, p, random_state=random_state)
        self.signal = signal
        # Generate model parameters
        self.beta_Z = self.signal * self.rng.standard_normal((self.p,self.K))

    def sample_X(self, n):
        X = self.rng.normal(0, 1, (n,self.p))
        return X.astype(np.float32)
    
    def compute_prob(self, X):
        f = np.matmul(X,self.beta_Z)
        prob = np.exp(f)
        prob_y = prob / np.expand_dims(np.sum(prob,1),1)
        return prob_y

    def sample_Y(self, X):
        prob_y = self.compute_prob(X)
        g = np.array([self.rng.multinomial(1,prob_y[i]) for i in range(X.shape[0])], dtype = float)
        classes_id = np.arange(self.K)
        y = np.array([np.dot(g[i],classes_id) for i in range(X.shape[0])], dtype = int)
        return y.astype(np.int32)

    
class DataModel_3(DataModel):
    def __init__(self, K, p, signal=1, random_state=None):
        super().__init__(K, p, random_state=random_state)
        self.signal = signal
        
    def sample_X(self, n):
        X = self.rng.normal(0, 1, (n,self.p))
        X[:,0] = self.rng.choice([-1,1], size=n, replace=True, p=[1/self.K, (self.K-1)/self.K])
        X[:,1] = self.rng.choice([-1,1], size=n, replace=True, p=[1/4, 3/4])
        X[:,2] = self.rng.choice([-1,1], size=n, replace=True, p=[1/2, 1/2])
        X[:,3] = self.rng.choice(self.K, n, replace=True)
        return X.astype(np.float32)
        
    def compute_prob(self, X):
        prob_y = np.zeros((X.shape[0], self.K))
        right_0 = X[:,0] > 0
        right_1 = X[:,1] > 0
        right_2 = X[:,2] > 0
        leaf_0 = np.where(1-right_0)[0]
        leaf_1 = np.where((right_0) * (1-right_1) * (1-right_2))[0]
        leaf_2 = np.where((right_0) * (1-right_1) * (right_2))[0]
        leaf_3 = np.where((right_0) * (right_1))[0]
        # Zeroth leaf: uniform distribution over all labels
        prob_y[leaf_0] = 1.0/self.K
        # First leaf: uniform distribution over the first half of the labels
        K_half = int(np.round(self.K/2))
        prob_y[leaf_1, 0:K_half] = 2.0/self.K
        prob_y[leaf_1, K_half:self.K] = 0
        # Second leaf: uniform distribution over the second half of the labels
        prob_y[leaf_2, 0:K_half] = 0
        prob_y[leaf_2, K_half:self.K] = 2.0/self.K
        # Third leaf: 90% probability to label determined by 4th variable
        X3 = np.round(X[leaf_3,3]).astype(int)
        for k in range(self.K):
            prob_y[leaf_3[X3==k],:] = (1-0.9)/(self.K-1.0)
            prob_y[leaf_3[X3==k],k] = 0.9
        # Standardize probabilities for each sample        
        prob_y = prob_y / prob_y.sum(axis=1)[:,None]
        return prob_y

    def sample_Y(self, X):
        prob_y = self.compute_prob(X)
        g = np.array([self.rng.multinomial(1,prob_y[i]) for i in range(X.shape[0])], dtype = float)
        classes_id = np.arange(self.K)
        y = np.array([np.dot(g[i],classes_id) for i in range(X.shape[0])], dtype = int)
        return y.astype(np.int32)
    

class DataModel_4(DataModel):
    def __init__(self, K, p, signal=1, imb=None, random_state=None):
        super().__init__(K, p, random_state=random_state)
        self.n_informative = np.maximum(1, int(self.p*0.5))
        self.class_sep = signal
        base = np.exp(-imb * np.arange(self.K))
        self.weights = base / np.sum(base)

    def sample(self, n):        
        X, Y = make_classification(n_samples=n, n_classes=self.K, weights=self.weights,
                                   n_features=self.p, n_informative=self.n_informative, 
                                   class_sep=self.class_sep, random_state=self.random_state)
        Y = Y.astype(np.int32)
        return X, Y

class DataModel_5(DataModel):
    # K-class classification where X|Y is heteroscedastic via a 2-component mixture:
    #  - "easy" component: tight Gaussian around class center
    #  - "hard" component: diffuse Gaussian, mean shifted toward other classes

    def __init__(
        self, K: int, p: int,
        pi_easy=0.6, delta_shift=0.40, center_scale=3.0, rho=None, epsilon0=0.01,
        # Parameters of the mixture
        sigma_easy=0.20, sigma_hard=1.20,
        random_state=None
    ):
        super().__init__(K, p, random_state=random_state)

        self.pi_easy = pi_easy
        self.sigma_easy = float(sigma_easy)
        self.sigma_hard = float(sigma_hard)
        self.delta_shift = float(delta_shift)
        self.center_scale = float(center_scale)
        self.epsilon0 = epsilon0

        if rho is None:
            self.rho = np.ones(self.K) / self.K
        else:
            rho = np.asarray(rho, dtype=float)
            self.rho = rho / rho.sum()

        # precompute class centers (mus) so sampling is stable/reproducible given seed
        self.mus = self._make_centers()

        # precompute "other-class means" for the boundary shift
        self.other_means = np.zeros_like(self.mus)
        for k in range(self.K):
            self.other_means[k] = self.mus[np.arange(self.K) != k].mean(axis=0)

        # store pi as vector
        if np.isscalar(self.pi_easy):
            self.pi_vec = np.full(self.K, float(self.pi_easy))
        else:
            pi = np.asarray(self.pi_easy, dtype=float)
            assert pi.shape == (self.K,)
            self.pi_vec = pi

    def _make_centers(self):
        #  - corners if K=4 and p>=2
        #  - otherwise random Gaussian centers scaled by center_scale

        if self.K == 4 and self.p >= 2:
            mus = np.zeros((self.K, self.p))
            c = self.center_scale
            mus[:, :2] = np.array([[+c, +c], [+c, -c], [-c, +c], [-c, -c]])
            if self.p > 2:
                mus[:, 2:] = self.rng.normal(0, 1.0, size=(self.K, self.p - 2))
            return mus
        else:
            return self.rng.normal(0, self.center_scale, size=(self.K, self.p))

    def sample_Y(self, n: int):
        # Sample labels with fixed class proportions rho
        return self.rng.choice(self.K, size=n, p=self.rho).astype(np.int32)

    def sample_X_given_Y(self, Y: np.ndarray):
        # Sample features conditional on given labels Y. Internally samples regime Z|Y (easy or hard) and then X|Y,Z.
        n = Y.shape[0]

        # regime: 1=easy, 0=hard
        Z = (self.rng.uniform(size=n) < self.pi_vec[Y]).astype(np.int32)

        mu_Y = self.mus[Y]
        mu_hard = mu_Y + self.delta_shift * (self.other_means[Y] - mu_Y)

        X = np.empty((n, self.p), dtype=float)

        idx_easy = (Z == 1)
        ne = idx_easy.sum()
        if ne > 0:
            X[idx_easy] = self.rng.normal(loc=mu_Y[idx_easy], scale=self.sigma_easy, size=(ne, self.p))

        idx_hard = ~idx_easy
        nh = idx_hard.sum()
        if nh > 0:
            X[idx_hard] = self.rng.normal(loc=mu_hard[idx_hard], scale=self.sigma_hard, size=(nh, self.p))

        return X

    def sample(self, n: int):
        Y0 = self.sample_Y(n)
        X = self.sample_X_given_Y(Y0)
        
        # Add intrinsic contamination
        T0 = contamination.construct_T_matrix_simple(self.K, self.epsilon0)
        contamination_process = contamination.LinearContaminationModel(T0, random_state=self.random_state)
        Y = contamination_process.sample_labels(Y0)     

        return X, Y
    
class DataModel_AP(DataModel):
    # K-class classification where X|Y is heteroscedastic via a 2-component mixture:
    #  - "easy" component: tight Gaussian around class center
    #  - "hard" component: diffuse Gaussian, mean shifted toward other classes

    def __init__(
        self, K: int, p: int,
        pi_easy=0.6, delta_shift=0.40, center_scale=3.0, rho=None, epsilon0=0.01,
        distribution_type="uniform",
        # Parameters of uniform
        length_interval_easy=None, length_interval_hard=None,
        # Parameters of truncated gaussian or "multivariate raised cosine"
        A_easy = None, R_easy = 1.0, A_hard = None, R_hard = 4.0,
        random_state=None):
        
        super().__init__(K, p, random_state=random_state)

        self.pi_easy = pi_easy
        self.delta_shift = delta_shift
        self.center_scale = center_scale
        self.epsilon0 = epsilon0

        if rho is None:
            self.rho = np.ones(self.K) / self.K
        else:
            rho = np.asarray(rho, dtype=float)
            self.rho = rho / rho.sum()

        self.distribution_type = distribution_type
        if distribution_type == "uniform":
            self.length_interval_easy = length_interval_easy
            self.length_interval_hard = length_interval_hard
        elif distribution_type == "truncated_gaussian":
            self.A_easy = A_easy
            self.R_easy = R_easy
            self.A_hard = A_hard
            self.R_hard = R_hard

        # precompute class centers (mus) so sampling is stable/reproducible given seed
        self.mus = self._make_centers()

        # precompute "other-class means" for the boundary shift
        self.other_means = np.zeros_like(self.mus)
        for k in range(self.K):
            self.other_means[k] = self.mus[np.arange(self.K) != k].mean(axis=0)

        # store pi as vector
        if np.isscalar(self.pi_easy):
            self.pi_vec = np.full(self.K, float(self.pi_easy))
        else:
            pi = np.asarray(self.pi_easy, dtype=float)
            assert pi.shape == (self.K,)
            self.pi_vec = pi

    def _make_centers(self):
        #  - corners if K=4 and p>=2
        #  - otherwise random Gaussian centers scaled by center_scale

        if self.K == 4 and self.p >= 2:
            mus = np.zeros((self.K, self.p))
            c = self.center_scale
            mus[:, :2] = np.array([[+c, +c], [+c, -c], [-c, +c], [-c, -c]])
            if self.p > 2:
                mus[:, 2:] = self.rng.normal(0, 1.0, size=(self.K, self.p - 2))
            return mus
        else:
            return self.rng.normal(0, self.center_scale, size=(self.K, self.p))

    def sample_Y(self, n: int):
        # Sample labels with fixed class proportions rho
        return self.rng.choice(self.K, size=n, p=self.rho).astype(np.int32)

    def sample_X_given_Y(self, Y: np.ndarray):
        # Sample features conditional on given labels Y. Internally samples regime Z|Y (easy or hard) and then X|Y,Z.
        n = Y.shape[0]

        Z = (self.rng.uniform(size=n) < self.pi_vec[Y]).astype(np.int32)

        mu_Y = self.mus[Y]
        mu_hard = mu_Y + self.delta_shift * (self.other_means[Y] - mu_Y)

        X = np.empty((n, self.p), dtype=float)

        idx_easy = (Z == 1)
        ne = idx_easy.sum()
        if ne > 0:
            if self.distribution_type=="uniform":
                X[idx_easy] = self.rng.uniform(low=mu_Y[idx_easy]-self.length_interval_easy/2, high=mu_Y[idx_easy]+self.length_interval_easy/2, size=(ne, self.p))
            elif self.distribution_type=="truncated_gaussian":
                X[idx_easy] = sample_truncated_gaussian(n=ne, R=self.R_easy, A=self.A_easy, mu=mu_Y[idx_easy])

        idx_hard = ~idx_easy
        nh = idx_hard.sum()
        if nh > 0:
            if self.distribution_type=="uniform":
                X[idx_hard] = self.rng.uniform(low=mu_hard[idx_hard]-self.length_interval_hard/2, high=mu_hard[idx_hard]+self.length_interval_hard/2, size=(nh, self.p))
            elif self.distribution_type=="truncated_gaussian":
                X[idx_hard] = sample_truncated_gaussian(n=nh, R=self.R_hard, A=self.A_hard, mu=mu_hard[idx_hard])
        return X

    def sample(self, n: int):
        Y0 = self.sample_Y(n)
        X = self.sample_X_given_Y(Y0)
        
        # Add intrinsic contamination
        T0 = contamination.construct_T_matrix_simple(self.K, self.epsilon0)
        contamination_process = contamination.LinearContaminationModel(T0, random_state=self.random_state)
        Y = contamination_process.sample_labels(Y0)     

        return X, Y


"""
class DataModel_draft(DataModel):
    def __init__(self, K, p, signal=1, random_state=None):
        super().__init__(K, p, random_state=random_state)
        self.signal = signal
        
    def sample_X(self, n):
        X = self.rng.normal(0, 1, (n,self.p))
        X[:,0] = self.rng.choice([-1,1], size=n, replace=True, p=[1/self.K, (self.K-1)/self.K])
        X[:,1] = self.rng.choice([-1,1], size=n, replace=True, p=[1/4, 3/4])
        X[:,2] = self.rng.choice([-1,1], size=n, replace=True, p=[1/2, 1/2])
        X[:,3] = self.rng.choice(self.K, n, replace=True)

        # Split observations into four leaves
        prob_y = np.zeros((X.shape[0], self.K))
        right_0 = X[:,0] > 0
        right_1 = X[:,1] > 0
        right_2 = X[:,2] > 0
        leaf_0 = np.where(1-right_0)[0]
        leaf_1 = np.where((right_0) * (1-right_1) * (1-right_2))[0]
        leaf_2 = np.where((right_0) * (1-right_1) * (right_2))[0]
        leaf_3 = np.where((right_0) * (right_1))[0]



        return X.astype(np.float32)
        
    def compute_prob(self, X):


        # Zeroth leaf: class 0 easily predictable, other classes have small separation
        # TO DO
        

        # First leaf: class 1 easily predictable, other classes have small separation
        # TO DO

        # Second leaf: class 2 easily predictable, other classes have small separation
        # TO DO

        # Second leaf: class 3 easily predictable, other classes have small separation
        # TO DO

        # Second leaf: uniform distribution over the second half of the labels
        prob_y[leaf_2, 0:K_half] = 0
        prob_y[leaf_2, K_half:self.K] = 2.0/self.K
        # Third leaf: 90% probability to label determined by 4th variable
        X3 = np.round(X[leaf_3,3]).astype(int)
        for k in range(self.K):
            prob_y[leaf_3[X3==k],:] = (1-0.9)/(self.K-1.0)
            prob_y[leaf_3[X3==k],k] = 0.9
        # Standardize probabilities for each sample        
        prob_y = prob_y / prob_y.sum(axis=1)[:,None]
        return prob_y

    def sample_Y(self, X):
        prob_y = self.compute_prob(X)
        g = np.array([self.rng.multinomial(1,prob_y[i]) for i in range(X.shape[0])], dtype = float)
        classes_id = np.arange(self.K)
        y = np.array([np.dot(g[i],classes_id) for i in range(X.shape[0])], dtype = int)
        return y.astype(np.int32)
"""

