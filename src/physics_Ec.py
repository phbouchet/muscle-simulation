"""
NumPy physics implementation of ACAP (as continuous as possible) energy, and its derivative
defined in 'EMU: Efficient Muscle Simulation in Deformation Space' by V. Modi et al.
"""

import numpy as np

def create_G(Q : np.ndarray, dim : tuple) -> np.ndarray:
    n, m = dim
    G = np.zeros((9*m, 3*n))

    for i in range (0,m):        
        Q1 = np.array([Q[i * 12], Q[i * 12 + 1], Q[i * 12 + 2]])
        Q2 = np.array([Q[i * 12 + 3], Q[i * 12 + 4] , Q[i * 12 + 5]])
        Q3 = np.array([Q[i * 12 + 6], Q[i * 12 + 7] , Q[i * 12 + 8]])
        Q4 = np.array([Q[i * 12 + 9], Q[i * 12 + 10], Q[i * 12 + 11]])

        R = np.linalg.inv((np.squeeze(np.array([Q1 - Q4, Q2 - Q4, Q3 - Q4])).T))
        R_11, R_21, R_31 = R[0,0], R[0,1], R[0,2]
        R_12, R_22, R_32 = R[1,0], R[1,1], R[1,2]
        R_13, R_23, R_33 = R[2,0], R[2,1], R[2,2]

        arr = np.array([[R_11,0,0,R_21,0,0,R_31,0,0,(-R_11-R_21-R_31),0,0],
                        [R_12,0,0,R_22,0,0,R_32,0,0,(-R_12-R_22-R_32),0,0],
                        [R_13,0,0,R_23,0,0,R_33,0,0,(-R_13-R_23-R_33),0,0]])

        G_i = np.tile(arr, (3,1))

        G_i[3:6] = np.roll(G_i[3:6], 1, 1)
        G_i[6:9] = np.roll(G_i[6:9], 2, 1)

        G[9*i:9*(i+1), 12*i:12*(i+1)] = G_i

    return G

def Ec(F : np.ndarray, G : np.ndarray, q : np.ndarray) -> int:
    return 0.5 * q.T @ G.T @ G @ q - q.T @ G.T @ F + 0.5 * F.T @ F

def argmin_Ec(F : np.ndarray, G : np.ndarray) -> np.ndarray:
    return np.linalg.inv(G.T @ G) @ G.T @ F


def d_Ec(F : np.ndarray, G : np.ndarray) -> np.ndarray:
    return -G @ Ec(F, G) + F