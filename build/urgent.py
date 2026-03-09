import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Beam momentum and target mass for GaN (adjust if needed)
p_beam = 3e9  # eV/c
M_target = 42 * 931.5e6  # GaN effective mass number × nucleon mass (eV)

# -------------------------------
# Combined Differential Cross-Section
# -------------------------------
files = {
    "pi-": "cross_section_pi.csv",
    "mu-": "cross_section_mu.csv",
    "e-":  "cross_section_e.csv"
}

plt.figure(figsize=(8,6))
for label, file in files.items():
    df_cs = pd.read_csv(file)
    plt.semilogy(df_cs["angle_deg"], df_cs["d_sigma_dOmega"],
                 drawstyle="steps-post", label=label)

plt.xlabel("Scattering angle (deg)")
plt.ylabel("dσ/dΩ (cm²/sr)")
plt.title("Differential Cross-Section — GaN")
plt.legend()
plt.grid(True, which="both", alpha=0.4)
plt.tight_layout()
plt.savefig("combined_cross_section_GaN.png", dpi=150)
plt.show()

# -------------------------------
# Combined RET vs Scattering Angle
# -------------------------------
files = {
    "pi-": "angular_correlation_pi.csv",
    "mu-": "angular_correlation_mu.csv",
    "e-":  "angular_correlation_e.csv"
}

plt.figure(figsize=(8,6))
for label, file in files.items():
    df_ac = pd.read_csv(file)
    plt.scatter(df_ac["Angle_deg"], df_ac["RET_eV"],
                s=8, alpha=0.6, label=label)

plt.xlabel("Scattering angle (deg)")
plt.ylabel("RET (eV)")
plt.title("RET vs Scattering Angle — GaN")
plt.yscale("log")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("combined_RET_vs_angle_GaN.png", dpi=150)
plt.show()

# -------------------------------
# Large-angle scattering fraction
# -------------------------------
particles = ["pi", "mu", "e"]
fractions = []
for p in particles:
    df_p = pd.read_csv(f"results_{p}.csv")
    frac = (df_p["Angle_deg"] > 5.0).sum() / len(df_p) * 100
    fractions.append(frac)

plt.figure(figsize=(6, 4))
plt.bar(particles, fractions, color=["steelblue", "orange", "green"])
plt.ylabel("Events with θ > 5° (%)")
plt.title("Large-angle scattering fraction by particle type (GaN, 3 GeV/c)")
plt.grid(axis="y", alpha=0.4)
plt.tight_layout()
plt.savefig("large_angle_fraction.png", dpi=150)
plt.show()

# -------------------------------
# RET distribution histogram
# -------------------------------
plt.figure(figsize=(7, 5))
for p, col in zip(particles, ["steelblue", "orange", "green"]):
    df_p = pd.read_csv(f"results_{p}.csv")
    theta = np.radians(df_p["Angle_deg"])
    RET = (2 * p_beam**2 / M_target) * (np.sin(theta/2)**2)
    plt.hist(RET[RET > 0], bins=50, log=True, alpha=0.5,
             label=p, color=col, density=True)

plt.axvline(x=20, color="red", linestyle="--", label="$E_d$ = 20 eV (GaN)")
plt.xlabel("Recoil energy transfer (eV)")
plt.ylabel("Normalised counts")
plt.title("RET distribution by particle type (GaN, 3 GeV/c)")
plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig("ret_distribution.png", dpi=150)
plt.show()


files = {
    "pi-": "angular_correlation_pi.csv",
    "mu-": "angular_correlation_mu.csv",
    "e-":  "angular_correlation_e.csv"
}

colors = {
    "pi-": "steelblue",
    "mu-": "orange",
    "e-":  "green"
}

plt.figure(figsize=(9,6))

for label, file in files.items():
    df = pd.read_csv(file)

    # Step plot for cumulative damage
    plt.step(df.index, df["Cumulative_Damage_eV"],
             where="post", lw=1.2, color=colors[label], label=f"{label} cumulative")

    # Scatter markers for large-angle events
    plt.scatter(df.index[df["LargeAngle"]],
                df["Cumulative_Damage_eV"][df["LargeAngle"]],
                s=12, color=colors[label], marker="o", edgecolor="red",
                alpha=0.8, zorder=3, label=f"{label} large-angle (θ > 5°)")

plt.xlabel("Event number")
plt.ylabel("Cumulative damage energy (eV)")
plt.title("Step-current style cumulative damage — GaN (π⁻, μ⁻, e⁻)")
plt.legend(fontsize=8, ncol=2)
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.savefig("combined_cumulative_damage_GaN.png", dpi=150)
plt.show()
