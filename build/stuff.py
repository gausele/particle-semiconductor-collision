import pandas as pd
import matplotlib.pyplot as plt

# Load your big CSV
df = pd.read_csv("summary.csv")

# --- Plot 1: Differential Cross-Section ---
plt.figure(figsize=(8,6))
for particle, group in df.groupby("particle"):
    plt.plot(group["angle_deg"], group["d_sigma_dOmega"],
             label=particle)
plt.xlabel("Scattering angle (deg)")
plt.ylabel("dσ/dΩ (cm²/sr)")
plt.yscale("log")
plt.title("Differential Cross-Section (combined)")
plt.legend()
plt.grid(True)
plt.savefig("combined_cross_section.png", dpi=300)
plt.show()

# --- Plot 2: RET vs Scattering Angle ---
plt.figure(figsize=(8,6))
for particle, group in df.groupby("particle"):
    plt.plot(group["Angle_deg"], group["avg_RET_eV"],
             label=particle)
plt.xlabel("Scattering angle (deg)")
plt.ylabel("Average RET (eV)")
plt.title("RET vs Scattering Angle (combined)")
plt.legend()
plt.grid(True)
plt.savefig("combined_RET_vs_angle.png", dpi=300)
plt.show()

