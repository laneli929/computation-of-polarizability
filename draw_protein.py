import matplotlib.pyplot as plt
from biopandas.pdb import PandasPdb

# 1. 读取 PDB 文件
# 这里的 'pdb/ca.pdb' 请替换为你的实际路径
ppdb = PandasPdb().read_pdb('pdb/ca.pdb')
atoms = ppdb.df['ATOM']

# 2. 提取坐标
x = atoms['x_coord'].values
y = atoms['y_coord'].values
z = atoms['z_coord'].values
elements = atoms['element_symbol'].values

# 3. 定义颜色映射函数
def get_color(el):
    mapping = {
        'C': [0.2, 0.2, 0.2], # 碳：深灰
        'N': [0.0, 0.0, 1.0], # 氮：蓝
        'O': [1.0, 0.0, 0.0], # 氧：红
        'S': [1.0, 1.0, 0.0]  # 硫：黄
    }
    return mapping.get(el, [0.0, 1.0, 0.0]) # 其他：绿

colors = [get_color(e) for e in elements]

# 4. 绘图
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 绘制散点图
scatter = ax.scatter(x, y, z, c=colors, s=20, alpha=0.8, edgecolors='none')

# 5. 设置坐标轴和标题
ax.set_xlabel('X (Å)')
ax.set_ylabel('Y (Å)')
ax.set_zlabel('Z (Å)')
ax.set_title('Protein colored by element')

# 保持比例一致 (等效于 MATLAB 的 axis equal)
# 注意：Matplotlib 3.6.0+ 支持 set_box_aspect('equal')
try:
    ax.set_box_aspect([1, 1, 1])
except AttributeError:
    pass

plt.show()