import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Load verified V2 dataset
DATA_PATH = Path("results/verified_v2_data.json")
with open(DATA_PATH, "r", encoding="utf-8") as f:
    data = json.load(f)

OUT = Path("poster_figures")
OUT.mkdir(exist_ok=True)
plt.rcParams.update({
    "font.size": 11.5,
    "axes.linewidth": 1.3,
    "font.sans-serif": "DejaVu Sans"
})

# ==============================================================================
# FIGURE 1: DME-Displacement Thermodynamics (V2 xTB/ALPB vs. TZVP/CPCM)
# ==============================================================================
swap = data["swap_energies"]

motifs = [
    "Alkoxy-OH\n(Methoxyethanol, CN=5)",
    "Amide C=O\n(DMA, CN=5)",
    "DME Reference\n([Mg(DME)₃]²⁺, CN=6)",
    "Zwitterion\n(Betaine, CN=6)"
]

xtb_vals = [
    swap["Methoxyethanol"]["xtb_alpb_kcal"],
    swap["DMA_amide"]["xtb_alpb_kcal"],
    0.00,
    swap["Betaine"]["xtb_alpb_kcal"]
]

dft_vals = [
    swap["Methoxyethanol"]["dft_tzvp_pcm_kcal"],
    swap["DMA_amide"]["dft_tzvp_pcm_kcal"],
    0.00,
    swap["Betaine"]["dft_tzvp_pcm_kcal"]
]

x = np.arange(len(motifs))
width = 0.36

fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=300)

ax.axhspan(-18, -3, color="#1b9e77", alpha=0.13, label='Optimal "Goldilocks" Desolvation Window')
ax.axhline(0, color="black", linewidth=1.1, linestyle="--")

rects1 = ax.bar(
    x - width / 2, xtb_vals, width,
    label="GFN2-xTB / ALPB(THF)",
    color="#7570b3", edgecolor="black", linewidth=1.1
)
rects2 = ax.bar(
    x + width / 2, dft_vals, width,
    label="B3LYP-D3BJ/def2-TZVP / CPCM(THF)",
    color="#1b9e77", edgecolor="black", linewidth=1.1
)

for rects in (rects1, rects2):
    for bar in rects:
        val = bar.get_height()
        if abs(val) < 1e-6:
            ax.text(bar.get_x() + bar.get_width() / 2, 0.8, "0.00",
                    ha="center", va="bottom", fontsize=9.5, color="#333333")
            continue
        y_pos = val + 0.8 if val > 0 else val - 0.9
        va = "bottom" if val > 0 else "top"
        weight = "bold" if val < 0 else "normal"
        ax.text(bar.get_x() + bar.get_width() / 2, y_pos, f"{val:+.2f}",
                ha="center", va=va, fontsize=10, fontweight=weight)

ax.annotate(
    "Sole Exothermic DME Stripper\n(Preserves Octahedral CN = 6)",
    xy=(3 + width / 2, dft_vals[3]),
    xytext=(1.55, -15.8),
    arrowprops=dict(arrowstyle="->", color="#0e5a43", lw=1.6),
    fontsize=9.8, fontweight="bold", color="#0e5a43",
    bbox=dict(boxstyle="round,pad=0.35", facecolor="#e8f8f3", edgecolor="#1b9e77")
)

ax.set_xticks(x)
ax.set_xticklabels(motifs, fontweight="bold", fontsize=10)
ax.set_ylabel(r"$\Delta E_{\mathrm{swap}}$ in Implicit THF (kcal/mol)", fontweight="bold")
ax.set_title("V2 Benchmark: Solvated DME-Displacement Thermodynamics", fontweight="bold", pad=12)
ax.set_ylim(-19.5, 28.5)
ax.grid(axis="y", linestyle=":", alpha=0.5)
ax.legend(loc="upper right", frameon=True, fontsize=9.8)

fig.tight_layout()
fig.savefig(OUT / "fig1_dme_displacement.png", dpi=300)
fig.savefig(OUT / "fig1_dme_displacement.svg")
plt.close(fig)

# ==============================================================================
# FIGURE 2: Relaxed Mg2+ Hopping Scans & Vertical Reduction Stability (V2)
# ==============================================================================
scans = data["relaxed_hopping_scans"]
red = data["vertical_reduction_alpb_thf"]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.7), dpi=300,
                               gridspec_kw={"width_ratios": [1.45, 1.0]})

# --- Panel A: 15-Point Relaxed Hopping Profiles ---
pb_fwd = scans["pure_betaine_fwd"]
hd_fwd = scans["hybrid_dma_fwd"]
pb_rev = scans["pure_betaine_rev"]
hd_rev = scans["hybrid_dma_rev"]

ax1.plot(
    pb_fwd["xs_valid"], pb_fwd["rel_e_kcal"],
    "o--", color="#d95f02", linewidth=2.1, markersize=5.5,
    label=f"Pure Betaine Fwd ($E_a = {pb_fwd['Ea_fwd_kcal']:.2f}$ kcal/mol)"
)
ax1.plot(
    hd_fwd["xs_valid"], hd_fwd["rel_e_kcal"],
    "s-", color="#1b9e77", linewidth=2.4, markersize=5.5,
    label=f"Betaine + DMA Bridge Fwd ($E_a = {hd_fwd['Ea_fwd_kcal']:.2f}$ kcal/mol)"
)
ax1.plot(
    pb_rev["xs_valid"], pb_rev["rel_e_kcal"],
    "^:", color="#e6ab02", linewidth=1.7, markersize=5.0, alpha=0.9,
    label=f"Pure Betaine Rev ($E_a = {pb_rev['Ea_fwd_kcal']:.2f}$ kcal/mol)"
)
ax1.plot(
    hd_rev["xs_valid"], hd_rev["rel_e_kcal"],
    "d-.", color="#2b5c8f", linewidth=1.8, markersize=5.0, alpha=0.9,
    label=f"Betaine + DMA Bridge Rev ($E_a = {hd_rev['Ea_fwd_kcal']:.2f}$ kcal/mol)"
)

delta_ea = hd_fwd["Ea_fwd_kcal"] - pb_fwd["Ea_fwd_kcal"]
pct_drop = abs(delta_ea) / pb_fwd["Ea_fwd_kcal"] * 100.0

ax1.annotate(
    "", xy=(0.0, hd_fwd["Ea_fwd_kcal"] + 0.5),
    xytext=(0.0, pb_fwd["Ea_fwd_kcal"] - 0.5),
    arrowprops=dict(arrowstyle="<->", color="#2b5c8f", lw=1.9)
)
ax1.text(
    0.16, 36.2,
    rf"$\Delta E_a = {delta_ea:.2f}\ \mathrm{{kcal/mol}}$" + f"\n({pct_drop:.1f}% Barrier Drop)",
    color="#1d436b", fontweight="bold", fontsize=9.5, va="center",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="#eef4fb", edgecolor="#2b5c8f")
)

ax1.set_xlabel(r"$\mathrm{Mg}^{2+}$ Reaction Coordinate ($\mathrm{\AA}$)", fontweight="bold")
ax1.set_ylabel("Relative Energy (kcal/mol)", fontweight="bold")
ax1.set_title("A. Relaxed 15-Point $\mathrm{Mg}^{2+}$ Hopping Scans", fontweight="bold", pad=10)
ax1.set_ylim(-3, 53)
ax1.grid(True, linestyle=":", alpha=0.5)
ax1.legend(loc="upper left", frameon=True, fontsize=8.4)

# --- Panel B: Vertical Electron Affinity (VEA) & Reduction Potential ---
species_keys = [
    "Control_[Mg(DME)3]2+",
    "Complex_[Mg(DME)2(Bet)]2+",
    "Complex_[Mg(DME)2(DMA)]2+"
]
species_labels = [
    "[Mg(DME)₃]²⁺\n(Control)",
    "[Mg(DME)₂(Bet)]²⁺\n(Zwitterion)",
    "[Mg(DME)₂(DMA)]²⁺\n(Amide)"
]
vea_vals = [red[k]["vea_ev"] for k in species_keys]
ered_vals = [red[k]["approx_E_red_vs_Mg_V"] for k in species_keys]
bar_colors = ["#7570b3", "#1b9e77", "#e6ab02"]

bars_b = ax2.bar(species_labels, vea_vals, width=0.52, color=bar_colors, edgecolor="black", linewidth=1.1)

for bar, vea, ered in zip(bars_b, vea_vals, ered_vals):
    ax2.text(
        bar.get_x() + bar.get_width() / 2, vea + 0.14,
        f"VEA: {vea:.2f} eV\n({ered:.2f} V vs Mg)",
        ha="center", va="bottom", fontsize=9.0, fontweight="bold"
    )

ax2.set_ylabel("Vertical Electron Affinity in ALPB(THF) (eV)", fontweight="bold")
ax2.set_title("B. Anodic Reduction Stability (ALPB-THF)", fontweight="bold", pad=10)
ax2.set_ylim(0, 7.1)
ax2.tick_params(axis="x", labelsize=9.2)
ax2.grid(axis="y", linestyle=":", alpha=0.5)

fig.tight_layout()
fig.savefig(OUT / "fig2_mg_migration_barrier.png", dpi=300)
fig.savefig(OUT / "fig2_mg_migration_barrier.svg")
plt.close(fig)

print("Generated updated V2 poster figures in poster_figures/:")
print(" - fig1_dme_displacement.png / .svg")
print(" - fig2_mg_migration_barrier.png / .svg")
