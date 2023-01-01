"""
NumPy physics implementation of Neo-Hookean fiber and isotropic functions,
as well as their derivatives and hessians.
"""

import numpy as np

def a(t : int) -> float:
    """ Linear activation function
    Args:
        t (int): Time variable
    """
    return (1/10) * t

def psi_fiber(F : np.ndarray, u : np.ndarray, t : int) -> np.ndarray:
    """ Eq.8
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        u (np.ndarray): Vertex-wise directional vector: (9m, 1)
        t (int): Time variable
    Returns:
        np.ndarray: (9m, 9m)
    """
    return a(t) * u @ F.T @ F @ u.T

def d_psi_fiber(F : np.ndarray, u : np.ndarray, t : int) -> np.ndarray:
    """ Derivative of Eq.8
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        u (np.ndarray): Vertex-wise directional vector: (9m, 1)
        t (int): Time variable
    Returns:
        np.ndarray: (9m, 9m)
    """
    return psi_fiber(F, u, t)

def psi_fiber_hessian(F : np.ndarray,  u : np.ndarray, t : int) -> np.ndarray:
    """ Hessian of Eq.8
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        u (np.ndarray): Vertex-wise directional vector: (9m, 1)
        t (int): Time variable
    Returns:
        np.ndarray: (9m, 9m)
    """
    return psi_fiber(F, u, t)

def psi_iso(F : np.ndarray):
    C = np.pi
    I = np.trace(F @ F.T)
    return C * (I - 3)

def d_psi_iso(F : np.ndarray):
    return 0

def psi_iso_hessian(F : np.ndarray):
    return 0
