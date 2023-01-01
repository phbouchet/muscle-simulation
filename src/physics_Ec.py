"""
NumPy physics implementation of ACAP (as continuous as possible) energy, and its derivative
defined in 'EMU: Efficient Muscle Simulation in Deformation Space' by V. Modi et al.
"""


import numpy as np

def E_c(F : np.ndarray, G : np.ndarray) -> np.ndarray:
    """ ACAP Energy Eq.4
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        G (np.ndarray): Deformation gradients according to q: (9m, 3n)
    Returns:
        np.ndarray: (3n, 1)
    """
    return (np.linalg.inv(G @ G.T) @ G).T @ F


def d_E_c(F : np.ndarray, G : np.ndarray) -> np.ndarray:
    """ Derivative of ACAP Energy Eq.10
    Args:
        F (np.ndarray): Deformation gradient: (9m, 1)
        G (np.ndarray): Deformation gradients according to q: (9m, 3n)
    Returns:
        np.ndarray: (9m, 1)
    """
    return G @ E_c(F, G) + F