import tidy3d as td
import numpy as np
import matplotlib.pyplot as plt

eps0 = 8.8541878128e-12  # 真空介电常数
c_const=3e8
lambda0=622e-9
pi=3.14
omega=c_const/lambda0*2*pi
eps_bj=1.0**2

# 1. 加载下载好的文件
sim_data = td.SimulationData.from_file('data/glucagon_h_2nm_air.hdf5')

# 1. 获取磁场分量 (确保监测器记录了 H)
Hx = sim_data["fieldmonitor_1"].Hx[:,:,:,0]
Hy = sim_data["fieldmonitor_1"].Hy[:,:,:,0]
Hz = sim_data["fieldmonitor_1"].Hz[:,:,:,0]

# 2. 计算旋度 J = curl(H)
# 注意：tidy3d 的数据是 xarray 格式，可以使用 .differentiate()
# Jx = dHz/dy - dHy/dz
Jx_total = Hz.differentiate("y") - Hy.differentiate("z")
# Jy = dHx/dz - dHz/dx
Jy_total = Hx.differentiate("z") - Hz.differentiate("x")
# Jz = dHy/dx - dHx/dy
Jz_total = Hy.differentiate("x") - Hx.differentiate("y")
# # 获取 E 的所有分量并插值 eps
# Ex = sim_data["fieldmonitor_1"].Ex[:,:,:,0]
# Ey = sim_data["fieldmonitor_1"].Ey[:,:,:,0]
# Ez = sim_data["fieldmonitor_1"].Ez[:,:,:,0]
#
# eps=sim_data["permittivitymonitor_2"].eps_xx.real[:,:,:,0]
# eps_aligned = eps.interp_like(Ex)
# print(eps_aligned)
# print(eps_aligned-1)
#
# Jx=-1j*omega*eps0*(eps_aligned-eps_bj)*Ex
# Jy=-1j*omega*eps0*(eps_aligned-eps_bj)*Ey
# Jz=-1j*omega*eps0*(eps_aligned-eps_bj)*Ez
Jx=Jx_total
Jy=Jy_total
Jz=Jz_total


x=Jx.x
y=Jy.y
z=Jz.z

# 计算磁矩密度分布 (Magnetization density)
mx_dens = 0.5 * (y * Jz - z * Jy)
my_dens = 0.5 * (z * Jx - x * Jz)
mz_dens = 0.5 * (x * Jy - y * Jx)

# 对 x, y, z 三个维度进行积分
m_x = mx_dens.integrate(("x", "y", "z")).values
m_y = my_dens.integrate(("x", "y", "z")).values
m_z = mz_dens.integrate(("x", "y", "z")).values

print(f"磁偶极矩分量 (Am^2):")
print(f"mx: {m_x}")
print(f"my: {m_y}")
print(f"mz: {m_z}")
#mx/Ex,things like that
alpha=np.abs(m_y)*1e-18/0.21

print('alpha_me=',alpha)