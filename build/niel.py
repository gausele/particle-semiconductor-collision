import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import mannwhitneyu

# ================================
# USER CONFIGURATION
# ================================
particle = "pi-"          # "pi+", "pi-", "mu+", "mu-", "e-"
material = "GaN"           # "Si", "GaAs", "GaN"

p_beam = 3e9              # beam momentum (eV/c)

masses = {
    "pi+": 139.6e6,
    "pi-": 139.6e6,
    "mu-": 105.7e6,
    "mu+": 105.7e6,
    "e-":  0.511e6
}

target_mass_numbers = {
    "Si":   28,
    "GaAs": 72,    # effective average (Ga=69.7, As=74.9)
    "GaN":  42     # effective average (Ga=69.7, N=14)
}

thresholds = {
    "Si":   25.0,  # eV
    "GaAs": 10.0,  # eV
    "GaN":  20.0   # eV
}

densities = {
    "Si":   2.33,  # g/cm³
    "GaAs": 5.32,
    "GaN":  6.15
}

thickness = 0.05   # cm
area     = 1.0     # cm²

# ================================
# LOAD GEANT4 RESULTS
# ================================
df = pd.read_csv("results_pi.csv")

if "Angle_deg" not in df.columns:
    raise ValueError("results.csv must contain column: Angle_deg")

# ================================
# MATERIAL + PARTICLE PROPERTIES
# ================================
m_particle = masses[particle]
M_target   = target_mass_numbers[material] * 931.5e6  # eV
E_th       = thresholds[material]                      # eV

density     = densities[material]
volume      = area * thickness                         # cm³
mass_target = density * volume                         # g

n_particles = len(df)
fluence     = n_particles / area                       # cm⁻²

# ================================
# RECOIL ENERGY TRANSFER (RET)
# Relativistic formula: T = (2p²/M) sin²(θ/2)
# Valid under M >> m (nuclear target >> projectile mass)
# ================================
theta = np.radians(df["Angle_deg"])
df["RET_eV"]   = (2 * p_beam**2 / M_target) * (np.sin(theta / 2)**2)
df["Displaced"] = df["RET_eV"] > E_th

n_displacements = df["Displaced"].sum()

print("=" * 40)
print("Simulation Summary")
print("=" * 40)
print(f"Particle  : {particle}")
print(f"Material  : {material}")
print(f"Particles : {n_particles}")
print(f"Threshold : {E_th} eV")
print(f"Displaced : {n_displacements}")
print(f"Avg RET   : {df['RET_eV'].mean():.4f} eV")

# ================================
# DAMAGE ENERGY PER EVENT
# Each displaced atom contributes its full RET;
# sub-threshold events contribute zero.
# ================================
df["DamageEnergy_eV"]    = df["RET_eV"] * df["Displaced"].astype(int)
df["Cumulative_Damage_eV"] = df["DamageEnergy_eV"].cumsum()

# Keep binary cumulative for reference
df["Cumulative_Displacements"] = df["Displaced"].cumsum()

# ================================
# NIEL ESTIMATION
# NIEL = total damage energy / (target mass × fluence)
# Units: MeV·cm²/g
# ================================
total_damage_energy = df["DamageEnergy_eV"].sum()
NIEL = (total_damage_energy / 1e6) / (mass_target * fluence)

print(f"NIEL      : {NIEL:.6e} MeV·cm²/g")

# ================================
# OBJECTIVE 1
# Differential Cross-Section dσ/dΩ
# dσ/dΩ = N(θ) / (Φ · ΔΩ)
# ================================
bins = np.linspace(0, 20, 40)
hist, edges = np.histogram(df["Angle_deg"], bins=bins)

solid_angle  = 2 * np.pi * (
    np.cos(np.radians(edges[:-1])) - np.cos(np.radians(edges[1:]))
)
d_sigma_dOmega = hist / (fluence * solid_angle)

plt.figure(figsize=(7, 5))
plt.semilogy(edges[:-1], d_sigma_dOmega, drawstyle="steps-post", color="steelblue")
plt.xlabel("Scattering angle (deg)")
plt.ylabel("dσ/dΩ (cm²/sr)")
plt.title(f"Differential cross-section — {material} ({particle})")
plt.grid(True, which="both", alpha=0.4)
plt.tight_layout()
plt.savefig("cross_section.png", dpi=150)
plt.show()

# ================================
# OBJECTIVE 2
# Experimental Damage Coefficient K_exp
# ΔI = K·Φ  →  K_exp = N_displacements / Φ
# ================================
K_exp = n_displacements / fluence
print(f"K_exp     : {K_exp:.6e} displacements·cm²")

# ================================
# OBJECTIVE 3
# Angular–Damage Correlation
# Plot cumulative DAMAGE ENERGY (not binary count)
# so large-angle events appear as visible step jumps
# ================================
df["LargeAngle"] = df["Angle_deg"] > 5.0

fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

# --- Top panel: cumulative damage energy ---
axes[0].plot(df.index, df["Cumulative_Damage_eV"],
             color="steelblue", lw=1, label="Cumulative damage energy")
axes[0].scatter(df.index[df["LargeAngle"]],
                df["Cumulative_Damage_eV"][df["LargeAngle"]],
                color="red", s=8, zorder=3,
                label="Large-angle events (θ > 5°)")
axes[0].set_ylabel("Cumulative damage energy (eV)")
axes[0].set_title(f"Angular–Damage Correlation — {material} ({particle})")
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.4)

# --- Bottom panel: per-event damage energy to show step sizes ---
axes[1].bar(df.index[df["Displaced"] & ~df["LargeAngle"]],
            df["DamageEnergy_eV"][df["Displaced"] & ~df["LargeAngle"]],
            width=1, color="steelblue", alpha=0.4, label="Small-angle damage")
axes[1].bar(df.index[df["Displaced"] & df["LargeAngle"]],
            df["DamageEnergy_eV"][df["Displaced"] & df["LargeAngle"]],
            width=1, color="red", alpha=0.7, label="Large-angle damage")
axes[1].set_xlabel("Event number")
axes[1].set_ylabel("Damage energy per event (eV)")
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.4)

plt.tight_layout()
plt.savefig("angular_damage_correlation.png", dpi=150)
plt.show()

# ================================
# STATISTICAL TEST
# Mann-Whitney U: tests whether large-angle events produce
# significantly higher damage energy than small-angle events.
# More physically meaningful than point-biserial on two binary variables.
# ================================
damage_small = df[~df["LargeAngle"]]["DamageEnergy_eV"]
damage_large = df[ df["LargeAngle"]]["DamageEnergy_eV"]

avg_small = damage_small.mean()
avg_large = damage_large.mean()
ratio     = avg_large / avg_small if avg_small > 0 else float("inf")

stat, pval = mannwhitneyu(damage_large, damage_small, alternative="greater")

print("=" * 40)
print("Objective 3 — Statistical Results")
print("=" * 40)
print(f"Avg damage energy, small-angle : {avg_small:.4f} eV")
print(f"Avg damage energy, large-angle : {avg_large:.4f} eV")
print(f"Ratio (large / small)          : {ratio:.1f}x")
print(f"Mann-Whitney U statistic       : {stat:.2f}")
print(f"p-value                        : {pval:.4e}")
if pval < 0.05:
    print("Result: Large-angle events produce SIGNIFICANTLY higher damage (p < 0.05)")
else:
    print("Result: No significant difference detected at p < 0.05")

large_angle_damage_fraction = damage_large.sum() / total_damage_energy \
                              if total_damage_energy > 0 else 0
large_angle_event_fraction  = df["LargeAngle"].sum() / n_particles

print(f"Large-angle events  : {100*large_angle_event_fraction:.2f}% of all events")
print(f"Large-angle damage  : {100*large_angle_damage_fraction:.2f}% of total damage energy")

# ================================
# SAVE RESULTS
# ================================
summary_data = {
    "particle":               [particle],
    "material":               [material],
    "n_particles":            [n_particles],
    "n_displacements":        [n_displacements],
    "avg_RET_eV":             [df["RET_eV"].mean()],
    "NIEL_MeV_cm2_g":        [NIEL],
    "K_exp":                  [K_exp],
    "avg_damage_small_eV":    [avg_small],
    "avg_damage_large_eV":    [avg_large],
    "damage_ratio":           [ratio],
    "MannWhitney_U":          [stat],
    "p_value":                [pval],
    "large_angle_event_frac": [large_angle_event_fraction],
    "large_angle_damage_frac":[large_angle_damage_fraction],
}

pd.DataFrame(summary_data).to_csv(
    "summary.csv",
    mode="a",
    header=not os.path.exists("summary.csv"),
    index=False
)

pd.DataFrame({
    "angle_deg":     edges[:-1],
    "d_sigma_dOmega": d_sigma_dOmega
}).to_csv("cross_section_pi.csv", index=False)

df[["Angle_deg", "RET_eV", "Displaced",
    "DamageEnergy_eV", "Cumulative_Damage_eV",
    "LargeAngle"]].to_csv("angular_correlation_pi.csv", index=False)

print("=" * 40)
print("Analysis complete. Files saved:")
print("  summary.csv")
print("  cross_section.csv")
print("  angular_correlation.csv")
print("  cross_section.png")
print("  angular_damage_correlation.png")
print("=" * 40)

plt.figure(figsize=(7,5))

plt.scatter(df["Angle_deg"], df["RET_eV"],
            s=8, alpha=0.6)

plt.xlabel("Scattering angle (deg)")
plt.ylabel("Recoil energy transferred (eV)")
plt.title(f"Damage vs Scattering Angle — {material} ({particle})")

plt.yscale("log")   # damage spans many orders of magnitude
plt.grid(True)

plt.show()

