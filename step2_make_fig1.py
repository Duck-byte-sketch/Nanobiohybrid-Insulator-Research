import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

v2 = json.loads(Path("results/verified_v2_data.json").read_text(encoding="utf-8"))
swap = v2.get("swap_energies", {})
fig_dir = Path("poster_figures")

order = (
    ("Methoxyethanol", "Methoxyethanol\n(Alcohol)"),
    ("DMA_amide", "DMA\n(Amide)"),
    ("Betaine", "Betaine\n(Zwitterion)"),
)
labels = tuple(pretty for k, pretty in order if k in swap)
xtb_vals = tuple(float(swap.get(k).get("xtb_alpb_kcal")) for k, _ in order if k in swap)
dft_vals = tuple(float(swap.get(k).get("dft_tzvp_pcm_kcal")) for k, _ in order if k in swap)
cns = tuple(int(swap.get(k).get("cn")) for k, _ in order if k in swap)

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=300)
x = np.arange(len(labels))
w = 0.34

b1 = ax.bar(x - w / 2, xtb_vals, width=w, label="GFN2-xTB (ALPB, THF)", color="#5b84b1", edgecolor="#222222", lw=0.7)
b2 = ax.bar(x + w / 2, dft_vals, width=w, label="DFT B3LYP-D3BJ / TZVP (PCM)", color="#1b4f72", edgecolor="#222222", lw=0.7)
ax.axhline(0.0, color="#222222", lw=1.0)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=11, fontweight="bold")
ax.set_ylabel("DME Displacement Energy (kcal/mol)", fontsize=11)

ymin, ymax = min(min(xtb_vals), min(dft_vals)) - 6.0, max(max(xtb_vals), max(dft_vals)) + 8.0
ax.set_ylim(ymin, ymax)

for i, (bar1, bar2, xv, dv, cn) in enumerate(zip(b1, b2, xtb_vals, dft_vals, cns)):
    ax.text(bar1.get_x() + w / 2, xv + (0.7 if xv >= 0 else -0.7), f"{xv:+.2f}", ha="center", va="bottom" if xv >= 0 else "top", fontsize=9.5)
    ax.text(bar2.get_x() + w / 2, dv + (0.7 if dv >= 0 else -0.7), f"{dv:+.2f}", ha="center", va="bottom" if dv >= 0 else "top", fontsize=9.5, fontweight="bold")
    badge_y = ymax - 2.4 if dv >= 0 else ymin + 1.5
    ax.text(x[i], badge_y, f"CN = {cn}", ha="center", va="center", fontsize=9.5, bbox=dict(boxstyle="round,pad=0.25", facecolor="#eef4f8", edgecolor="#8aa8be", lw=0.8))

ax.legend(loc="upper right", frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=10)
ax.set_title("Figure 1: Thermodynamics of DME Displacement from [Mg(DME)3]2+", fontsize=13, fontweight="bold", pad=18)

fig.tight_layout()
fig.savefig(fig_dir / "fig1_dme_displacement.png", dpi=300, bbox_inches="tight", pad_inches=0.3)
fig.savefig(fig_dir / "fig1_dme_displacement.svg", bbox_inches="tight", pad_inches=0.3)
plt.close(fig)
print("Step 2 complete: Saved unclipped Figure 1.")
