# computing alpha_ee or alpha_em
# set the external field
# use the same direction for the dipole moment

import tidy3d as td
import numpy as np

eps0 = 8.8541878128e-12  # Vacuum permittivity

# 1. Load the downloaded simulation file
sim_data = td.SimulationData.from_file('5nmsphere_demo.hdf5')

# Get all E-field components
Ex = sim_data["fieldmonitor_1"].Ex[:,:,:,0]
Ey = sim_data["fieldmonitor_1"].Ey[:,:,:,0]
Ez = sim_data["fieldmonitor_1"].Ez[:,:,:,0]

# eps=sim_data["permittivitymonitor_2"].eps_xx.real[:,:,:,0]
# eps_aligned = eps.interp_like(Ex)
#
# Dx = Ex * eps_aligned
# Dy = Ey * eps_aligned
# Dz = Ez * eps_aligned

divE = Ex.differentiate("x") + Ey.differentiate("y") + Ez.differentiate("z")
#
# total_charge = eps0 * divE.integrate(coord=['x', 'y', 'z'])
# print(total_charge)

# Get the coordinates of divE
x_coords = divE.x
y_coords = divE.y
z_coords = divE.z

# print(x_coords)
# test

# Calculate the three components of the dipole moment (px, py, pz)
# We need to multiply the coordinates by their corresponding divergence values
px_density = x_coords * divE
py_density = y_coords * divE
pz_density = z_coords * divE

px = eps0 * px_density.integrate(coord=['x', 'y', 'z'])
py = eps0 * py_density.integrate(coord=['x', 'y', 'z'])
pz = eps0 * pz_density.integrate(coord=['x', 'y', 'z'])

print(f"Dipole Moment: p_x={px.values}, p_y={py.values}, p_z={pz.values}")

# use px or py? be careful!
E_ext=74 # using external field, 74 in this case
alpha = np.abs(px.values) * 1e-18 / E_ext

print('numerical result of alpha:', alpha)

## 4. Analytical solution for polarizability alpha (sphere)
pi = 3.14
eps_r = 1.6**2
eps_medium = 1.33**2
R = 5e-9 # Theoretical radius
alpha_th = 4 * pi * eps0 * R**3 * (eps_r - eps_medium) / (eps_r + 2 * eps_medium)
print("analytical result of alpha:", alpha_th)