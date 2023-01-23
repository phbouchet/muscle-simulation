import os
import numpy as np

from mesh_tools import *
from physics_lib import *

def EMU(F : np.ndarray, Q : np.ndarray, u : np.ndarray, dim : tuple, mesh_path : str, iter : int = None):
    dt = 0.5
    t, i = 0, 0
    σ = 10
    c, p, α = 1e-4, 0.5, 10
    G = create_G(Q, dim)
    e_i = E(F, G, u, t, α)
    g = d_E(F, G, u, t, α)

    epsilon_1 = -1e-4
    epsilon_2 = -1e-3
    
    while ((np.linalg.norm(g) > epsilon_1) or (E(F, G, u, t, α) - e_i) < epsilon_2):

        if (iter != None and i > iter):
            break

        e_i = E(F, G, u, t, α)
        g   = d_E(F, G, u, t, α)
        H_inv = E_hessian_inv(F, G, u, t, α, dim)
        d = H_inv @ g
        
        E_ = e_i

        iteration = 0
        F_temp = np.copy(F)

        while (E_ < e_i + σ * c * g.T @ d):
            F_temp = F + σ * d
            q = argmin_Ec(F_temp, G)
            σ *= p

            E_ = E(F_temp, G, u, t, α)
            
            out_path = f"{os.getcwd()}/data/EMU_mesh_{i}_{iteration}.obj"
            export_mesh(q, mesh_path, out_path)

            if (iteration > 100000):
                print("=== Diverged, halting execution ===")
                return
            
            iteration += 1

        print(f"Iteration: {i} | Converged at: {iteration}")
        F = F_temp + σ * d
        i += 1
        t += dt