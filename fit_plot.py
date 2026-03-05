import numpy as np
import matplotlib.pyplot as plt

# 1. Setup Data
# Base materials used for the fit
materials = ["Glucagon", "Calmodulin", "Rubisco","myo"]
mw_base = np.array([3.49, 16.85, 279.81,17.87])
kappa_base = np.array([0.010, 0.0016, 0.00039,0.008])

# # New data point (to be compared against the fit)
# myo_name = "Myoglobin"
# myo_mw = 17.87
# myo_kappa = 0.008

# 2. Perform Log-Log Fitting (κ = A * MW^n)
# log10(kappa) = n * log10(MW) + log10(A)
log_mw = np.log10(mw_base)
log_kappa = np.log10(kappa_base)

coeffs = np.polyfit(log_mw, log_kappa, 1)
n, logA = coeffs
A = 10**logA

# 3. Generate Fit Curve for plotting
mw_fit_range = np.linspace(min(mw_base)*0.8, max(mw_base)*1.1, 500)
kappa_fit_curve = A * (mw_fit_range**n)

# 4. Visualization
plt.figure(figsize=(8, 6))

# Plot Fit Line
fit_label = f'Fit (Base): κ = {A:.2e} · MW^{{ {n:.2f} }}'
plt.plot(mw_fit_range, kappa_fit_curve, '--', color='tab:blue', alpha=0.7, label=fit_label)

# Plot Base Data Points
plt.scatter(mw_base, kappa_base, color='tab:red', marker='s', s=80, label='Base Materials', zorder=5)
#
# # Plot Myoglobin
# plt.scatter(myo_mw, myo_kappa, color='green', marker='o', s=100, label=f'{myo_name} (Test)', zorder=5)

# Annotations for Base Materials
for name, x, y in zip(materials, mw_base, kappa_base):
    plt.text(x, y * 1.15, f"{name}\n({x}, {y:.1e})", ha='center', fontsize=9)

# # Annotation for Myoglobin
# plt.text(myo_mw, myo_kappa * 1.15, f"{myo_name}\n({myo_mw}, {myo_kappa:.1e})",
#          ha='center', fontsize=9, color='green', fontweight='bold')

# Formatting
plt.xlabel("Molecular Weight (kDa)")
plt.ylabel("Chiral Parameter κ")
plt.title("Protein Chiral Parameter: Trend vs. Observation")
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.yscale('linear') # Keeping linear scale as per your second request

plt.tight_layout()
plt.show()

# 5. Output Results
# expected_kappa = A * (myo_mw**n)
print(f"Fit Equation: κ = {A:.4e} * MW^({n:.4f})")
# print(f"Myoglobin MW: {myo_mw}")
# print(f"Expected κ (from fit): {expected_kappa:.4e}")
# print(f"Observed κ: {myo_kappa:.4e}")
# print(f"Deviation: {((myo_kappa - expected_kappa) / expected_kappa) * 100:.2f}%")