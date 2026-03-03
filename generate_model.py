import numpy as np
from Bio.PDB import PDBParser
import trimesh
import os  # 新增：用于处理文件路径
# -----------------------------
# 1. 读取 PDB 文件
# -----------------------------
pdb_file = "pdb/ca.pdb"  # 改成你的 PDB 文件路径
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)
base_name = os.path.splitext(os.path.basename(pdb_file))[0]
atom_coords = []

for atom in structure.get_atoms():
    # 获取原子坐标
    pos = atom.get_coord()  # 返回 [x, y, z] 的 numpy array
    atom_coords.append(pos)

atom_coords = np.array(atom_coords)
print(f"读取了 {len(atom_coords)} 个原子坐标")

# -----------------------------
# 2. 生成小球 STL
# -----------------------------
sphere_radius = 2  # 半径 1 Å -> 直径 2 Å
sphere_mesh = trimesh.creation.icosphere(subdivisions=2, radius=sphere_radius)

all_meshes = []

for coord in atom_coords:
    # 平移小球到原子坐标
    mesh_copy = sphere_mesh.copy()
    mesh_copy.apply_translation(coord)
    all_meshes.append(mesh_copy)

# 合并所有小球
combined = trimesh.util.concatenate(all_meshes)

# # -----------------------------
# # 3. 输出 STL 文件
# # -----------------------------
# output_file = "protein_atoms.stl"
# combined.export(output_file)
# print(f"STL 文件已生成: {output_file}")

# -----------------------------
# 3. 输出 STL 文件 (动态命名)
# -----------------------------
# 将 base_name 加在原定文件名前面
output_file = f"{base_name}_protein_atoms.stl"

combined.export(output_file)
print(f"STL 文件已生成: {output_file}")