import shutil
import numpy as np

def import_mesh(mesh_path : str):
    vertex_data = []
    face_data = []
    for line in open(mesh_path, "r"):
        if line.startswith('#'):
            continue
        values = line.split()
        if not values:
            continue
        if values[0] == 'v':
            v = list(map(float, values[1:4]))
            vertex_data.append(v)
        elif values[0] == 'f':
            f = list(map(lambda x: int(x.split('/')[0]),  values[1:4]))
            face_data.append(f)
    vertices = np.array(vertex_data)
    faces = np.array(face_data)
    return vertices, faces 

def export_mesh(mesh : np.ndarray, base_file : str, out_file : str):
    shutil.copy(base_file, out_file)

    with open(base_file, 'r') as file:
        data = file.readlines()

    tet_i = 0
    l = len(data)
    for i in range(l):
        if data[i][0] != 'v' or data[i][1] != ' ':
            continue
        data[i] = f"v {mesh[tet_i*3][0]} {mesh[tet_i*3+1][0]} {mesh[tet_i*3+2][0]}\n"
        tet_i += 1

    with open(out_file, 'w') as out_file:
        out_file.writelines(data)