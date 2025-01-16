import numpy as np
from scipy.stats import norm
from sklearn.model_selection import train_test_split
from sklearn.datasets import make_classification
import pdb

def sigmoid(x):
    return(1/(1 + np.exp(-x)))


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
        X, Y = make_classification(n_samples=n, n_classes=self.K, n_features=self.p, n_informative=self.n_informative, 
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
    def __init__(self, K, p, imbalance_ratios=None, imb=None, random_state=None):
        """
        Parameters:
        - K: Number of classes.            
        - p: Number of features.
        - imbalance_ratios: Array of length K with relative proportions of each class.
        - random_state: Seed for reproducibility.
        """
        super().__init__(K, p, random_state=random_state)

        if imbalance_ratios is not None:
            self.imbalance_ratios = np.array(imbalance_ratios)/np.sum(imbalance_ratios)
        elif imb is not None:
            base = np.exp(-imb * np.arange(K))
            self.imbalance_ratios = base / np.sum(base)
        else:
            raise ValueError("Either imbalance_ratios or imb must be provided.")
        
    def sample(self, n):        
        X, Y = make_classification(n_samples=n, n_classes=self.K, n_features=self.p, n_informative=self.n_informative, 
                                   class_sep=self.class_sep, random_state=self.random_state)
        Y = Y.astype(np.int32)
        return X, Y

    def sample_X(self, n):
        return self.rng.normal(0, 1, size=(n, self.p))

    def sample_Y(self, X):
        """Sample Y with the specified class imbalance."""
        return self.rng.choice(self.K, size=X.shape[0], p=self.imbalance_ratios)

class DataModel_5(DataModel):
    def __init__(self, K, p, signal=1, imbalance_ratios=None, imb=None, random_state=None):
        super().__init__(K, p, random_state=random_state)
        self.n_informative = np.maximum(1, int(self.p*0.5))
        self.class_sep = signal

        if imbalance_ratios is not None:
            self.imbalance_ratios = np.array(imbalance_ratios)/np.sum(imbalance_ratios)
        elif imb is not None:
            base = np.exp(-imb * np.arange(K))
            self.imbalance_ratios = base / np.sum(base)
        else:
            raise ValueError("Either imbalance_ratios or imb must be provided.")

    def sample(self, n):        
        X, Y = make_classification(
            n_samples=n, 
            n_classes=self.K, 
            n_features=self.p, 
            n_informative=self.n_informative, 
            class_sep=self.class_sep, 
            random_state=self.random_state
        )
        Y = Y.astype(np.int32)
    
        # Post-sampling adjustment to control class balance
        adjusted_indices = []
        for k in range(self.K):
            # Indices of samples belonging to class k
            class_k_indices = np.where(Y == k)[0]
        
            # Number of samples to retain for class k
            num_samples = int(self.imbalance_ratios[k] * n)
        
            # Select a subset of indices for the desired class balance
            if num_samples <= len(class_k_indices):
                adjusted_indices.extend(self.rng.choice(class_k_indices, size=num_samples, replace=False))
            else:
                # Oversample and add slight noise to avoid duplicates
                oversampled_indices = self.rng.choice(class_k_indices, size=num_samples - len(class_k_indices), replace=True)
                adjusted_indices.extend(class_k_indices)  # Include all available samples
                adjusted_indices.extend(oversampled_indices)

                 # Optional: Add noise to duplicate samples to make them unique
                oversampled_indices = np.array(oversampled_indices)
                noise = self.rng.normal(scale=0.01, size=X[oversampled_indices].shape)  # Small Gaussian noise
                X[oversampled_indices] += noise
    
        adjusted_indices = np.array(adjusted_indices)
        return X[adjusted_indices], Y[adjusted_indices]

