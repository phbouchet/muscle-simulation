"""
NumPy physics implementation the numeric solution of the scalar energy minimization problem in Eq.1
defined in 'EMU: Efficient Muscle Simulation in Deformation Space' by V. Modi et al.
"""

import numpy as np
from physics_psi import *
from physics_Ec  import *

def E(F : np.ndarray, G : np.ndarray, u : np.ndarray, t : int, α : int) -> np.ndarray:
    """ Discretized energy Eq.5
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        G (np.ndarray): Deformation gradients according to q: (9m, 3n)
        u (np.ndarray): Vertex-wise directional vector: (9m, 1)
        t (int): Time variable
        α (int): Hyperparameter
    Returns:
        np.ndarray: (9m, 9m)
    """
    return psi_iso(F) + psi_fiber(F, u, t) + α * E_c(F, G)

def d_E(F : np.ndarray, G : np.ndarray, u : np.ndarray, t : int, α : int) -> np.ndarray:
    """ Derivative of discretized energy Eq.9
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        G (np.ndarray): Deformation gradients according to q: (9m, 3n)
        u (np.ndarray): Vertex-wise directional vector: (9m, 1)
        t (int): Time variable
        α (int): Hyperparameter
    Returns:
        np.ndarray: (9m, 9m)
    """
    return d_psi_iso(F) + d_psi_fiber(F, u, t) + α * d_E_c(F, G)

def E_hessian(F : np.ndarray, G : np.ndarray, u : np.ndarray, t : int, α : int, dim : tuple) -> np.ndarray:
    """ Hessian of discretized energy Eq.9
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        G (np.ndarray): Deformation gradients according to q: (9m, 3n)
        u (np.ndarray): Vertex-wise directional vector: (9m, 1)
        t (int): Time variable
        α (int): Hyperparameter
        dim (tuple): Tuple containing values of n and m
    Returns:
        np.ndarray: (9m, 9m)
    """
    n, m = dim
    
    k = 3*n

    w, v = np.linalg.eig(G.T @ G)
    eig_val, eig_vec = np.array([w]), np.squeeze(np.array([v]))

    if (3*n > 48):
        k = 48 # Justified after Eq.13
        eig_val = np.sort(eig_val)[:k] # k smallest eigenvalues

    Φ = eig_vec # (k, 3n) or (3n, 3n)
    Λ = np.multiply(np.identity(k), eig_val) # (k, k)
    
    I = np.identity(9*m)
    H = psi_iso_hessian(F) + psi_fiber_hessian(F, u, t) + α * I  #(9m, 9m)

    B = Φ @ G.T

    # Hessian matrix of E
    hessian = H - α * B.T @ np.linalg.inv(Λ) @ B

    return hessian

def E_hessian_inv(F : np.ndarray, G : np.ndarray, u : np.ndarray, t : int, α : int, dim : tuple) -> np.ndarray:
    """ Inverse hessian of discretized energy Eq.9
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        G (np.ndarray): Deformation gradients according to q: (9m, 3n)
        u (np.ndarray): Vertex-wise directional vector: (9m, 1)
        t (int): Time variable
        α (int): Hyperparameter
        dim (tuple): Tuple containing values of n and m
    Returns:
        np.ndarray: (9m, 9m)
    """
    n, m = dim 
    k = 3*n

    w, v = np.linalg.eig(G.T @ G)
    eig_val, eig_vec = np.array([w]), np.squeeze(np.array([v]))

    if (3*n > 48):
        k = 48 # Justified after Eq.13
        eig_val = np.sort(eig_val)[:k] # k smallest eigenvalues

    Φ = eig_vec # (k, 3n) or (3n, 3n)
    Λ = np.multiply(np.identity(k), eig_val) # (k, k)
    
    I = np.identity(9*m)
    H = psi_iso_hessian(F) + psi_fiber_hessian(F, u, t) + α * I #(9m, 9m)
    H_inv = np.linalg.inv(H) # (9m, 9m)

    B = Φ @ G.T # (k, 9m)
    
    # Inverse of hessian expression
    hessian_inv = H_inv + α * H_inv @ B.T @ Λ - α*B @ H_inv @ B.T @ B @ H_inv

    return hessian_inv