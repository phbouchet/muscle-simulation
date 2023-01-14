import os
import numpy as np

from mesh_tools import *
from physics_lib import *

def EMU(F : np.ndarray, Q : np.ndarray, u : np.ndarray, dim : tuple, mesh_path : str):
    t = 0
    sigma = 10
    c, p, α = 1e-4, 0.5, 100
    G = create_G(Q, dim)
    g = d_E(F, G, u, t, α)

    epsilon_1 = -1e-4
    epsilon_2 = -1e-3
    
    while (np.linalg.norm(g) > epsilon_1 and t < 10):
        F_temp = np.copy(F)
        e_i = E(F, G, u, t, α)
        g = d_E(F, G, u, t, α)
        H_inv = E_hessian_inv(F, G, u, t, α, dim)
        d = H_inv @ g
        
        E_ = e_i

        iteration_conv = 0

        while (E_ < e_i + sigma * c * g.T @ d):
            F_temp = F + sigma * d
            sigma = p * sigma
            E_ = E(F_temp, G, u, t, α)
            q = argmin_Ec(F_temp, G)
            export_mesh(q, mesh_path, f"{os.getcwd()}/data/EMU_mesh_{t}_{iteration_conv}.obj")

            if (iteration_conv > 100000):
                print("=== Diverged, halting execution ===")
                return
            
            iteration_conv += 1


        print(f"Iteration: {t} | Converged at: {iteration_conv}")
        F = F_temp + sigma * d
        t += 1