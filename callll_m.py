import tidy3d as td
import numpy as np

# --- 1. Constants and Parameters ---
eps0 = 8.8541878128e-12
c_const = 299792458
lambda0 = 622e-9
omega = 2.0 * np.pi * c_const / lambda0
E_ext = 74  # External field amplitude

# --- 2. Load Data ---
sim_data = td.SimulationData.from_file('data/ca_march.hdf5')

# Get E-field components (assuming monitor recorded at index 0 for frequency)
Ex = sim_data["fieldmonitor_1"].Ex.isel(f=0)
Ey = sim_data["fieldmonitor_1"].Ey.isel(f=0)
Ez = sim_data["fieldmonitor_1"].Ez.isel(f=0)

# Get Permittivity and align with E-field grid
eps = sim_data["permittivitymonitor_2"].eps_xx.real.isel(f=0)
eps_aligned = eps.interp_like(Ex)

# --- 3. Calculate Induced Current Density J ---
# J = -i * omega * P = -i * omega * (eps_r - eps_bg) * eps0 * E
# Assuming background permittivity eps_bg = 1.0
eps_bg = 1.0
chi = (eps_aligned )

Jx = -1j * omega * eps0 * chi * Ex
Jy = -1j * omega * eps0 * chi * Ey
Jz = -1j * omega * eps0 * chi * Ez

# --- 4. Define Coordinates ---
# Tidy3D DataArrays have coordinates accessible via .x, .y, .z
x = Jx.x
y = Jy.y
z = Jz.z

# --- 5. Calculate Magnetic Dipole Moment Density ---
# m = 1/2 * integral( r x J ) dV
# mx = 0.5 * (y*Jz - z*Jy)
# my = 0.5 * (z*Jx - x*Jz)
# mz = 0.5 * (x*Jy - y*Jx)

mx_density = 0.5 * (y * Jz - z * Jy)
my_density = 0.5 * (z * Jx - x * Jz)
mz_density = 0.5 * (x * Jy - y * Jx)

# --- 6. Integrate to find m ---
m_x = mx_density.integrate(coord=['x', 'y', 'z'])
m_y = my_density.integrate(coord=['x', 'y', 'z'])
m_z = mz_density.integrate(coord=['x', 'y', 'z'])

print(f"Magnetic Dipole Moment (Am^2):")
print(f"mx: {m_x.values}")
print(f"my: {m_y.values}")
print(f"mz: {m_z.values}")

# --- 7. Calculate Magnetic Polarizability alpha_mm ---
# Magnetizability is often defined as m = alpha_mm * H_ext
# Or sometimes m = alpha * E_ext/c. Here we follow your specific scaling:
alpha_mm = np.abs(m_x.values) * 1e-18 / E_ext

print(f"Numerical result of alpha_mm: {alpha_mm}")