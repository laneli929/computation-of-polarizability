#create a model according to PDB file
#warning : not watertight !!!!!!

import numpy as np
from Bio.PDB import PDBParser
from skimage import measure
import trimesh
from scipy.ndimage import distance_transform_edt
import os

PDB_FILE = "pdb/1gcn.pdb"
STL_FILE = os.path.splitext(os.path.basename(PDB_FILE))[0] + "_fast_dense.stl"

# PDB_FILE = "pdb/ca.pdb"
# STL_FILE = "ca_fast_dense.stl"

VOXEL_SIZE = 1.0
ATOM_RADIUS = 2.0
PADDING = 5.0

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", PDB_FILE)

coords = []
for atom in structure.get_atoms():
    if atom.element != "H":
        coords.append(atom.coord)
coords = np.array(coords)

print(f"Loaded {len(coords)} atoms")

min_corner = coords.min(axis=0) - ATOM_RADIUS - PADDING
max_corner = coords.max(axis=0) + ATOM_RADIUS + PADDING

xs = np.arange(min_corner[0], max_corner[0], VOXEL_SIZE)
ys = np.arange(min_corner[1], max_corner[1], VOXEL_SIZE)
zs = np.arange(min_corner[2], max_corner[2], VOXEL_SIZE)

nx, ny, nz = len(xs), len(ys), len(zs)
print(f"Voxel grid: {nx} × {ny} × {nz}")

# 1. Create an empty voxel grid (False = empty)
grid = np.zeros((nx, ny, nz), dtype=bool)

# 2. Map atom coordinates to the nearest voxel points (single point occupation)
idx = np.floor((coords - min_corner) / VOXEL_SIZE).astype(int)
idx = np.clip(idx, 0, np.array([nx-1, ny-1, nz-1]))

grid[idx[:,0], idx[:,1], idx[:,2]] = True

# 3. Calculate distance field (fast)
dist = distance_transform_edt(~grid) * VOXEL_SIZE

# 4. Marching cubes algorithm
verts, faces, normals, _ = measure.marching_cubes(
    dist,
    level=ATOM_RADIUS,
    spacing=(VOXEL_SIZE, VOXEL_SIZE, VOXEL_SIZE)
)
verts += min_corner

print(f"Generated mesh: {len(verts)} vertices, {len(faces)} faces")

# 5. Clean mesh
mesh = trimesh.Trimesh(verts, faces, process=False)
# mesh.remove_duplicate_faces()
# mesh.remove_degenerate_faces()
# mesh.remove_unreferenced_vertices()
# mesh.fill_holes()
# mesh.fix_normals()

print("Watertight:", mesh.is_watertight)

mesh.export(STL_FILE)
print(f"Saved clean STL: {STL_FILE}")