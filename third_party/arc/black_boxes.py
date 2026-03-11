import numpy as np
from sklearn import svm
from sklearn import ensemble
from sklearn import calibration
from sklearn.neural_network import MLPClassifier

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

import copy


class Oracle:
    def __init__(self, model):
        self.model = model
    
    def fit(self,X,y):
        return self

    def predict(self, X):
        return self.model.sample(X)        

    def predict_proba(self, X):
        if(len(X.shape)==1):
            X = X.reshape((1,X.shape[0]))
        prob = self.model.compute_prob(X)
        prob = np.clip(prob, 1e-6, 1.0)
        prob = prob / prob.sum(axis=1)[:,None]
        return prob


class SVC:
    def __init__(self, calibrate=False,
                 kernel = 'linear',
                 C = 1,
                 clip_proba_factor = 0.1,
                 random_state = 2020):
        self.model = svm.SVC(kernel = kernel,
                             C = C,
                             probability = True,
                             random_state = random_state)
        self.calibrate = calibrate
        self.num_classes = 0
        self.factor = clip_proba_factor
        
    def fit(self, X, y):
        self.num_classes = len(np.unique(y)) 
        self.model_fit = self.model.fit(X, y)
        if self.calibrate:
            self.calibrated = calibration.CalibratedClassifierCV(self.model_fit,
                                                                 method='sigmoid',
                                                                 cv=10)
        else:
            self.calibrated = None
        return copy.deepcopy(self)

    def predict(self, X):
        return self.model_fit.predict(X)

    def predict_proba(self, X):        
        if(len(X.shape)==1):
            X = X.reshape((1,X.shape[0]))
        if self.calibrated is None:
            prob = self.model_fit.predict_proba(X)
        else:
            prob = self.calibrated.predict_proba(X)
        prob = np.clip(prob, self.factor/self.num_classes, 1.0)
        prob = prob / prob.sum(axis=1)[:,None]
        return prob

class RFC:
    def __init__(self, calibrate=False,
                 n_estimators = 1000,
                 criterion="gini", 
                 max_depth=None,
                 max_features=None,
                 min_samples_leaf=1,
                 clip_proba_factor=0.1,
                 random_state = 2020):
        
        self.model = ensemble.RandomForestClassifier(n_estimators=n_estimators,
                                                     criterion=criterion,
                                                     max_depth=max_depth,
                                                     max_features=max_features,
                                                     min_samples_leaf=min_samples_leaf,
                                                     random_state = random_state)
        self.calibrate = calibrate
        self.num_classes = 0
        self.factor = clip_proba_factor
        
    def fit(self, X, y):
        self.num_classes = len(np.unique(y)) 
        self.model_fit = self.model.fit(X, y)
        if self.calibrate:
            self.calibrated = calibration.CalibratedClassifierCV(self.model_fit,
                                                                 method='sigmoid',
                                                                 cv=10)
        else:
            self.calibrated = None
        return copy.deepcopy(self)

    def predict(self, X):
        return self.model_fit.predict(X)

    def predict_proba(self, X):        
        if(len(X.shape)==1):
            X = X.reshape((1,X.shape[0]))
        if self.calibrated is None:
            prob = self.model_fit.predict_proba(X)
        else:
            prob = self.calibrated.predict_proba(X)
        prob = np.clip(prob, self.factor/self.num_classes, 1.0)
        prob = prob / prob.sum(axis=1)[:,None]
        return prob

class NNet:
    def __init__(self, calibrate=False,
                 hidden_layer_sizes = 64,
                 batch_size = 128,
                 learning_rate_init = 0.01,
                 max_iter = 20,
                 clip_proba_factor = 0.1,
                 random_state = 2020):
        
        self.model = MLPClassifier(hidden_layer_sizes=hidden_layer_sizes,
                                   batch_size=batch_size,
                                   learning_rate_init=learning_rate_init,
                                   max_iter=max_iter,
                                   random_state=random_state)
        self.calibrate = calibrate
        self.num_classes = 0
        self.factor = clip_proba_factor
        
    def fit(self, X, y):
        self.num_classes = len(np.unique(y)) 
        self.model_fit = self.model.fit(X, y)
        if self.calibrate:
            self.calibrated = calibration.CalibratedClassifierCV(self.model_fit,
                                                                 method='sigmoid',
                                                                 cv=10)
        else:
            self.calibrated = None
        return copy.deepcopy(self)

    def predict(self, X):
        return self.model_fit.predict(X)

    def predict_proba(self, X):        
        if(len(X.shape)==1):
            X = X.reshape((1,X.shape[0]))
        if self.calibrated is None:
            prob = self.model_fit.predict_proba(X)
        else:
            prob = self.calibrated.predict_proba(X)
        prob = np.clip(prob, self.factor/self.num_classes, 1.0)
        prob = prob / prob.sum(axis=1)[:,None]
        return prob

class MLPBlackBox:
    def __init__(self,
                 hidden_sizes=(256, 128),
                 dropout=0.3,
                 lr=1e-3,
                 batch_size=64,
                 max_epochs=100,
                 patience=10,          # early stopping
                 clip_proba_factor=1e-5,
                 random_state=2020):

        self.hidden_sizes = hidden_sizes
        self.dropout = dropout
        self.lr = lr
        self.batch_size = batch_size
        self.max_epochs = max_epochs
        self.patience = patience
        self.factor = clip_proba_factor
        self.random_state = random_state
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.num_classes = 0
        self.model_fit = None

    def _build_model(self, input_dim, num_classes):
        layers = []
        in_dim = input_dim
        for h in self.hidden_sizes:
            layers += [nn.Linear(in_dim, h), nn.BatchNorm1d(h), nn.ReLU(), nn.Dropout(self.dropout)]
            in_dim = h
        layers.append(nn.Linear(in_dim, num_classes))
        return nn.Sequential(*layers)

    def fit(self, X, y):
        torch.manual_seed(self.random_state)
        self.num_classes = len(np.unique(y))
        input_dim = X.shape[1]

        model = self._build_model(input_dim, self.num_classes).to(self.device)

        # Split off 10% for early stopping validation
        n_val = max(1, int(0.1 * len(X)))
        idx = np.random.default_rng(self.random_state).permutation(len(X))
        X_val, y_val = X[idx[:n_val]], y[idx[:n_val]]
        X_tr,  y_tr  = X[idx[n_val:]], y[idx[n_val:]]

        # Convert to tensors
        def to_tensors(Xa, ya):
            return TensorDataset(torch.tensor(Xa, dtype=torch.float32),
                                 torch.tensor(ya, dtype=torch.long))

        train_loader = DataLoader(to_tensors(X_tr, y_tr),
                                  batch_size=self.batch_size, shuffle=True)
        X_val_t = torch.tensor(X_val, dtype=torch.float32).to(self.device)
        y_val_t = torch.tensor(y_val, dtype=torch.long).to(self.device)

        optimizer = optim.Adam(model.parameters(), lr=self.lr, weight_decay=1e-4)
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5)
        criterion = nn.CrossEntropyLoss()

        best_val_loss = float('inf')
        best_state = None
        patience_counter = 0

        for epoch in range(self.max_epochs):
            # Train
            model.train()
            for Xb, yb in train_loader:
                Xb, yb = Xb.to(self.device), yb.to(self.device)
                optimizer.zero_grad()
                loss = criterion(model(Xb), yb)
                loss.backward()
                optimizer.step()

            # Validate
            model.eval()
            with torch.no_grad():
                val_loss = criterion(model(X_val_t), y_val_t).item()
            scheduler.step(val_loss)

            if val_loss < best_val_loss:
                best_val_loss = val_loss
                best_state = copy.deepcopy(model.state_dict())
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= self.patience:
                    break

        model.load_state_dict(best_state)
        model.eval()
        self.model_fit = model
        return copy.deepcopy(self)

    def predict_proba(self, X):
        if len(X.shape) == 1:
            X = X.reshape(1, -1)
        X_t = torch.tensor(X, dtype=torch.float32).to(self.device)
        with torch.no_grad():
            logits = self.model_fit(X_t)
            prob = torch.softmax(logits, dim=1).cpu().numpy()
        prob = np.clip(prob, self.factor / self.num_classes, 1.0)
        prob = prob / prob.sum(axis=1, keepdims=True)
        return prob

    def predict(self, X):
        return np.argmax(self.predict_proba(X), axis=1)


