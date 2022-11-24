import os
import trimesh

import numpy as np
import scipy as sci

from mesh_tools import *

print(f"Loading file from : {os.getcwd()}/data/arm_tet.obj")

mesh_path = f"{os.getcwd()}/data/arm_tet.obj"

vertices, faces = load_obj_mesh(mesh_path)
mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
triangle = get_triangle(vertices, faces[0])
print(triangle)