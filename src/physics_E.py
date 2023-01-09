"""
NumPy physics implementation the numeric solution of the scalar energy minimization problem in Eq.1
defined in 'EMU: Efficient Muscle Simulation in Deformation Space' by V. Modi et al.
"""

import numpy as np
from physics_Ec  import *
from physics_psi import *

def E(F : np.ndarray, G : np.ndarray, u : np.ndarray, t : int, α : int) -> int:
    q = argmin_Ec(F, G)
    return psi_iso(F) + psi_fiber(F, u, t) + α * Ec(F, G, q)[0,0]

def d_E(F : np.ndarray, G : np.ndarray, u : np.ndarray, t : int, α : int) -> np.ndarray:
    return d_psi_iso(F) + d_psi_fiber(F, u, t) + α * d_Ec(F, G)

def E_hessian(F : np.ndarray, G : np.ndarray, u : np.ndarray, t : int, α : int, dim : tuple) -> np.ndarray:
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
    n, m = dim 
    k = 3*n

    w, v = np.linalg.eig(np.linalg.inv(G.T @ G))
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
    hessian_inv = H_inv + α * H_inv @ B.T @ np.linalg.inv(Λ - α*B @ H_inv @ B.T) @ B @ H_inv

    return hessian_inv