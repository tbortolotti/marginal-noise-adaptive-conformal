import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
import pdb
import cvxpy as cp

def estimate_c_const(n, n_mc=1000):
    R = np.zeros((n_mc,))
    for b in range(n_mc):
        U = np.sort(np.random.uniform(0,1,size=(n,)))
        R[b] = np.max(np.arange(1,n+1)/n - U)
    c = np.mean(R)
    return c

def eval_delta_marg_a(beta_0, betas, W, n):
    K = W.shape[0]
    # Identity matrix of size K
    I_K = np.eye(K)

    # Constructing a matrix with the i-th row as beta_i
    Beta_matrix = cp.vstack([betas[i] * cp.Constant(np.ones(K)/K) for i in range(K)])

    # Constructing Omega matrix
    Omega = W - beta_0 * I_K - Beta_matrix

    # Constants for objective function
    const_a = np.sqrt(np.log(K*n+1))

    # Calculate c(n)
    c_n = estimate_c_const(n)

    # Loss functions
    loss_1 = c_n * (beta_0 + cp.sum(betas)/K)
    loss_2_a = cp.norm(Omega, 1) * const_a
    
    loss = loss_1 + 2/np.sqrt(n) * loss_2_a
    return loss

def eval_delta_marg_b(beta_0, betas, W, n):
    K = W.shape[0]

    # Identity matrix of size K
    I_K = np.eye(K)

    # Constructing a matrix with the i-th row as beta_i
    Beta_matrix = cp.vstack([betas[i] * cp.Constant(np.ones(K)/K) for i in range(K)])

    # Constructing Omega matrix
    Omega = W - beta_0 * I_K - Beta_matrix

    # Constants for objective function
    const_b = 24 * ((2 * np.log(K)+1)/(2 * np.log(K)-1)) * np.sqrt(2*K*np.log(K))

    # Calculate c(n)
    c_n = estimate_c_const(n)

    # Loss functions
    loss_1 = c_n * (beta_0 + cp.sum(betas)/K)
    loss_2_b = cp.norm(Omega, "inf") * const_b
    
    loss = loss_1 + 2/np.sqrt(n) * loss_2_b
    return loss

def eval_delta_marg(beta_0, betas, W, n):
    loss_a = eval_delta_marg_a(beta_0, betas, W, n).value
    loss_b = eval_delta_marg_b(beta_0, betas, W, n).value
    return np.minimum(loss_a, loss_b)

def solve_optim_problem_a(W, n):
    K = W.shape[0]
    # Variables
    beta_0 = cp.Variable()
    betas = cp.Variable(K)
   
    # Solve problem A
    loss = eval_delta_marg_a(beta_0, betas, W, n)
    objective = cp.Minimize(loss)
    problem = cp.Problem(objective, [])
    problem.solve()

    # Return results
    optim_value = problem.value
    optim_beta_0 = beta_0.value
    optim_betas = betas.value
    return optim_value
    
def solve_optim_problem_b(W, n):
    K = W.shape[0]
    # Variables
    beta_0 = cp.Variable()
    betas = cp.Variable(K)

    # Solve problem A
    loss = eval_delta_marg_b(beta_0, betas, W, n)
    objective = cp.Minimize(loss)
    problem = cp.Problem(objective, [])
    problem.solve()
    
    # Return results
    optim_value = problem.value
    optim_beta_0 = beta_0.value
    optim_betas = betas.value
    return optim_value

def eval_delta_marg_opt(W, n):
    value_a = solve_optim_problem_a(W, n)
    value_b = solve_optim_problem_b(W, n)
    
    if value_a <= value_b:
        value = value_a
    else:
        value = value_b
        
    return value
