import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import pdb
from functools import lru_cache
import sys

from cln.utils_ecdf import EmpiricalCDF2DGrid

def construct_grid(num_steps, grid_type='uniform'):

    if grid_type == 'uniform':
        grid = np.linspace(0, 1, num_steps)
    elif grid_type == 'centered':
        # Step 1: Create a linear grid
        linear_grid = np.linspace(0, 1, np.ceil(num_steps/2).astype(int))
        # Step 2: Apply a quadratic transformation
        # This transformation makes the intervals smaller around 0.5
        grid_1 = 0.5 * (np.cos(np.pi * (linear_grid - 0.5)))
        grid_2 = 1-0.5 * (np.cos(np.pi * (linear_grid - 0.5)))
        grid = np.concatenate([grid_1,grid_2])
        grid = np.sort(grid)
    else:
        raise ValueError("Unsupported grid type. Choose 'uniform' or 'centered'.")

    return grid

def compute_joint_F_hat(scores, t_values):
    joint_F_hat_vals = dict()
    K = max(key[0] for key in scores.keys()) + 1

    for t1 in t_values:
        for t2 in t_values:
            if (t1, t2) not in joint_F_hat_vals:
                joint_F_hat_vals[(t1, t2)] = dict()
            for l in range(K):
                for k in range(K):
                    s_k_l = scores[(l,k)]
                    for k_prime in range(K):
                        s_k_prime_l = scores[(l,k_prime)]
                        joint_F_hat_vals[(t1, t2)][(l, k, k_prime)] = np.mean((s_k_l <= t1) & (s_k_prime_l <= t2))

    return joint_F_hat_vals

def cov_empirical_process(grid, W, scores, F_hat, n_cal):
    K = max(key[0] for key in scores.keys()) + 1
    n_grid = len(grid)

    Sigma = np.zeros((n_grid, n_grid))
    exp1 = np.zeros((n_grid,))

    for l in range(K):
        n_l = len(scores[(l,0)])
        for k in range(K):
            exp1 += W[k,l] * n_l/n_cal * F_hat[(l,k)](grid)
            for k_prime in range(K):
                # Joint ecdf
                scores_x = scores[(l,k)]
                scores_y = scores[(l,k_prime)]
                ecdf_2d_grid = EmpiricalCDF2DGrid(scores_x, scores_y)
                joint_F_hat = ecdf_2d_grid.evaluate_grid(grid, grid)
                # Update sum
                Sigma += W[k,l] * W[k_prime,l] * n_l/n_cal * joint_F_hat

    Sigma = Sigma - np.dot(exp1.reshape(-1,1), exp1.reshape(-1,1).T)

    return Sigma


def simulate_gaussian_process(grid, num_samples, n_cal, scores, W, F_hat):
    # Compute the covariance matrix
    n_grid = len(grid)

    Sigma = cov_empirical_process(grid, W, scores, F_hat, n_cal)

    # Add a small value to the diagonal for numerical stability
    Sigma += 1e-6 * np.eye(n_grid)

    # Mean vector
    mu = np.zeros((n_grid,))

    # Generate samples
    samples = np.random.multivariate_normal(mu, Sigma, size=(num_samples))

    return samples

def simulate_supremum(h, num_samples, grid_type, n_cal, scores, W, F_hat):
    num_steps = int(np.ceil(1.0 / h))
    grid = construct_grid(num_steps, grid_type)

    samples = simulate_gaussian_process(grid, num_samples, n_cal, scores, W, F_hat)
    suprema = np.max(samples,1)

    expected_supremum = np.mean(suprema)
    return expected_supremum

@lru_cache(maxsize=None)
def richardson_recursive(h, func, j, k):
    """Recursive function for Richardson extrapolation."""
    if j == 0:
        return func(h)
    if k == 0:
        return func(h/np.power(2,j))
    else:
        R1 = richardson_recursive(h, func, j, k-1)
        R2 = richardson_recursive(h, func, j-1, k-1)
        return (np.power(2.0,k) * R1 - R2) / (np.power(2.0,k) - 1.0)

def richardson_extrapolation(h, func, max_j):
    """Perform Richardson extrapolation using recursion and memoization."""
    return richardson_recursive(h, func, max_j, max_j)
