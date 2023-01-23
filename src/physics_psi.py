"""
NumPy physics implementation of Neo-Hookean fiber and isotropic functions,
as well as their derivatives and hessians.
"""

import numpy as np

def a(t : int) -> float:
    return 1/1000 * t

def psi_fiber(F : np.ndarray, u : np.ndarray, t : int) -> int:
    return a(t) * (F.T @ F @ u.T @ u)[0,0]

def d_psi_fiber(F : np.ndarray, u : np.ndarray, t : int) -> np.ndarray:
    return a(t) * (F @ u.T @ u)

def psi_fiber_hessian(F : np.ndarray,  u : np.ndarray, t : int, m : int) -> np.ndarray:
    A = np.zeros((9*m, 9*m))
    
    for i in range(0,m):
        u_i = u[3*i:3*(i+1)] @ u[3*i:3*(i+1)].T
        A[9*i:9*(i+1), 9*i:9*(i+1)] = np.tile(u_i, (3,3))

    return a(t) * A

def psi_iso(F : np.ndarray):
    return 0

def d_psi_iso(F : np.ndarray):
    return 0

def psi_iso_hessian(F : np.ndarray):
    return 0
