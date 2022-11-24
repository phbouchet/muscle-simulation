import numpy as np

def load_obj_mesh(mesh_path):
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

def get_triangle(vertices : np.ndarray, faces : np.ndarray) -> np.ndarray:
    triangle = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            triangle[i, j] = vertices[faces[i] - 1, j]
    return triangle