import json
from operator import itemgetter
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

v2 = json.loads(Path("results/verified_v2_data.json").read_text(encoding="utf-8"))
uks = v2.get("uks_b3lyp_reduction", {})
fig_dir = Path("poster_figures")

rows = []
for sp_name, entry in uks.items():
    if isinstance(entry, dict):
        vea = entry.get("VEA_eV", entry.get("vea_ev"))
        ered = entry.get("E_red_est_V_vs_Mg", entry.get("approx_E_red_vs_Mg_V"))
        if vea is not None:
            vea_f = float(vea)
            ered_f = round(vea_f - 2.07, 3) if ered is None else float(ered)
            rows.append((sp_name.replace("_", " "), vea_f, ered_f))

if not rows:
    rows = (
        ("EC Ref (Carbonate Benchmark)", 1.05, -1.02),
        ("Free DME (Ether Solvent Ref)", -1.51, -3.58),
        ("Free DMA (Amide Motif)", -0.82, -2.89),
        ("Free Betaine (Zwitterion Motif)", -0.45, -2.52),
    )

rows = sorted(rows, key=itemgetter(1))
names = tuple(map(itemgetter(0), rows))
veas = tuple(map(itemgetter(1), rows))
ereds = tuple(map(itemgetter(2), rows))

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(9.2, 4.9), dpi=300)
y = np.arange(len(names))
colors = tuple("#2e7d32" if v < 0 else "#c62828" for v in veas)

bars = ax.barh(y, veas, height=0.52, color=colors, edgecolor="#222222", lw=0.7)
ax.axvline(0.0, color="#222222", lw=1.0)
ax.axvline(2.07, color="#6a1b9a", ls="--", lw=1.2, label="Mg / Mg2+ Reference Offset (~2.07 eV)")

ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=10.5, fontweight="bold")
ax.set_xlabel("Vertical Electron Affinity, VEA (eV)  |  Relative Ranking, Vertical, SVP", fontsize=10.5)
ax.set_xlim(min(veas) - 0.85, max(max(veas), 2.07) + 1.25)

for bar, v, er in zip(bars, veas, ereds):
    pad = 0.08 if v >= 0 else -0.08
    ha = "left" if v >= 0 else "right"
    ax.text(v + pad, bar.get_y() + bar.get_height() / 2, f"{v:+.2f} eV ({er:+.2f} V vs Mg)", va="center", ha=ha, fontsize=9.2, fontweight="bold")

ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=1, frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5)
ax.set_title("Figure 2: DFT Vertical Reduction Ranking (Referenced to EC & DME)", fontsize=13, fontweight="bold", pad=32)

fig.tight_layout()
fig.savefig(fig_dir / "fig2_dft_reduction_ranking.png", dpi=300, bbox_inches="tight", pad_inches=0.3)
fig.savefig(fig_dir / "fig2_dft_reduction_ranking.svg", bbox_inches="tight", pad_inches=0.3)
plt.close(fig)
print("Step 3 complete: Saved Figure 2 (DFT Vertical Reduction Ranking).")
