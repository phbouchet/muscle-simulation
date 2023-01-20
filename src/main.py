#!/usr/bin/env python3

import os
import numpy as np
import argparse

from mesh_tools import *
from physics_EMU import *

np.random.seed(42)

def main(mesh_file : str, iterations : int = None):
    mesh_path = f"{os.getcwd()}/{mesh_file}"
    vertices, _ = import_mesh(mesh_path)

    n = len(vertices) 
    m = int(n/4)
    dim = (n, m)

    F = np.tile(np.ravel(np.identity(3))[:, None], (m, 1))
    Q = np.ravel(np.array(vertices))[:, None]
    u = np.ones((3*m, 1))

    EMU(F, Q, u, dim, mesh_path, iterations)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Input parameters for EMU algorithm.")
    parser.add_argument('file', type=str, help="Target mesh file in obj format for simulation.")
    parser.add_argument('-i', '--iterations', type=int, help="Number of iterations")
    
    args = parser.parse_args()

    mesh_file = args.file
    iterations = args.iterations

    main(mesh_file, iterations)