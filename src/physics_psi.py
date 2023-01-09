"""
NumPy physics implementation of Neo-Hookean fiber and isotropic functions,
as well as their derivatives and hessians.
"""

import numpy as np

def a(t : int) -> float:
    return (1/10) * t

def psi_fiber(F : np.ndarray, u : np.ndarray, t : int) -> int:
    return a(t) * (F.T @ F @ u.T @ u)[0,0]

def d_psi_fiber(F : np.ndarray, u : np.ndarray, t : int) -> np.ndarray:
    return a(t) * (F @ u.T @ u)

def psi_fiber_hessian(F : np.ndarray,  u : np.ndarray, t : int) -> np.ndarray:
    return psi_fiber(F, u, t)

def psi_iso(F : np.ndarray):
    return 0

def d_psi_iso(F : np.ndarray):
    return 0

def psi_iso_hessian(F : np.ndarray):
    return 0
