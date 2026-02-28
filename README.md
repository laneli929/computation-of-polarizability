# Polarizability Computation Suite

This tool provides an efficient numerical method for calculating the **electric ($\alpha_{ee}$)** and **magnetic ($\alpha_{em}$)** polarizability of biomolecules or nanoparticles using Tidy3D simulation data.

---

## 1. Simulation Setup
* **Monitor Configuration:** Ensure `fieldmonitor_1` (or your specific monitor name) is set to capture the full 3D field distribution around the structure.
* **Excitation Strategy:** Use two symmetric plane waves to isolate specific polarizabilities. Adjust the source settings to maintain only the required components (e.g., $E_x$ or $H_y$).

## 2. Extraction & Analysis
* **Dipole Calculation:** The script computes the induced dipole moments ($p_x, p_y, p_z$) by integrating the E-field divergence across the simulation volume.
* **Polarizability Derivation:** * **$\alpha_{ee}$**: Calculated as the ratio $p_x / E_{ext}$.
    * **$\alpha_{em}$**: Calculated as the ratio $p_y / H_{ext}$.
* **Field Verification:** The external field values ($E_{ext}, H_{ext}$) can be extracted from the simulation result.
## 3. Validation & Benchmarking
To ensure numerical accuracy, the script's output is compared against the **Clausius-Mossotti** analytical solution for a dielectric sphere ($R = 5\text{nm}$, $n=1.6$) in a medium ($n=1.33$).

| Method | $\alpha_{ee}$ Value ($\text{F}\cdot\text{m}^2$) |
| :--- | :--- |
| **Analytical (Theory)** | $1.803 \times 10^{-36}$ |
| **Numerical (Tidy3D)** | $1.797 \times 10^{-36}$ |
| **Relative Error** | **~0.33%** |

---
