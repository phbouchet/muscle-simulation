import os
import numpy as np

from EMU import *
from mesh_tools import *

np.random.seed(42)

mesh_file = "simple_tet.obj"
mesh_path = f"{os.getcwd()}/data/{mesh_file}"
vertices, faces = load_mesh(mesh_path)

n = len(vertices)
m = int(n / 4)
dim = (n, m)

F = np.ravel(np.identity(3))[:, None]
Q = np.ravel(np.array(vertices))[:, None]
u = np.ones((3*m, 1))

EMU(F, Q, u, dim, mesh_path)