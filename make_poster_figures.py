from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

OUT = Path("poster_figures")
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 12, 'axes.linewidth': 1.3})

# --- Figure 1: DME Displacement ---
fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=300)
motifs = ['Alkoxy-OH\n(Cellulose)', 'Amide C=O\n(Silk/Peptide)', 'DME\n(Reference)', 'Zwitterion\n(Betaine)', 'Carboxylate\n(Alginate)']
energies = [23.41, 11.69, 0.00, -12.34, -31.14]
colors = ['#d95f02', '#e6ab02', '#7570b3', '#1b9e77', '#e7298a']

ax.axhspan(-22, -4, color='#1b9e77', alpha=0.14, label='Optimal "Goldilocks" Desolvation Window')
ax.axhline(0, color='black', linewidth=1.0, linestyle='--')
bars = ax.bar(motifs, energies, color=colors, edgecolor='black', linewidth=1.1, width=0.55)

for bar, val in zip(bars, energies):
    y_pos = val + 1.5 if val >= 0 else val - 3.2
    weight = 'bold' if abs(val + 12.34) < 0.1 else 'normal'
    ax.text(bar.get_x() + bar.get_width()/2, y_pos, f'{val:+.2f}',
            ha='center', va='bottom' if val >= 0 else 'top', fontsize=11, fontweight=weight)

ax.annotate('DFT (B3LYP-D3BJ/def2-SVP):\n-36.92 kcal/mol (Gas Phase)',
            xy=(3, -12.34), xytext=(1.45, -27.0),
            arrowprops=dict(arrowstyle='->', color='#1b9e77', lw=1.6),
            fontsize=10, fontweight='bold', color='#0e5a43',
            bbox=dict(boxstyle='round,pad=0.35', facecolor='#e8f8f3', edgecolor='#1b9e77'))

ax.set_ylabel(r'$\Delta E_{\mathrm{swap}}$ in ALPB(THF) (kcal/mol)', fontweight='bold')
ax.set_title('DME-Displacement Thermodynamics Across Bio-Hybrid Motifs', fontweight='bold', pad=12)
ax.set_ylim(-38, 31)
ax.grid(axis='y', linestyle=':', alpha=0.5)
ax.legend(loc='upper right', frameon=True, fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "fig1_dme_displacement.png", dpi=300)
fig.savefig(OUT / "fig1_dme_displacement.svg")
plt.close(fig)

# --- Figure 2: Mg2+ Migration Barrier ---
fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=300)
xs = np.array([-1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75])
e_pure = np.array([0.00, 17.23, 34.90, 45.81, 46.38, 36.50, 19.76, 3.23])
e_hyb  = np.array([0.00,  7.05, 17.06, 24.33, 25.43, 20.15, 11.11, 3.89])

ax.plot(xs, e_pure, 'o--', color='#d95f02', linewidth=2.2, markersize=7,
        label='Pure Zwitterion Channel ($E_a = 46.38$ kcal/mol)')
ax.plot(xs, e_hyb, 's-', color='#1b9e77', linewidth=2.5, markersize=7,
        label='Betaine + Amide C=O Bridge ($E_a = 25.43$ kcal/mol)')

ax.annotate('', xy=(0.25, 25.6), xytext=(0.25, 46.1),
            arrowprops=dict(arrowstyle='<->', color='#2b5c8f', lw=2.0))
ax.text(0.35, 35.8, r'$\Delta E_a = -20.95\ \mathrm{kcal/mol}$' + '\n(45% Barrier Reduction)',
        color='#1d436b', fontweight='bold', fontsize=10.5, va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#eef4fb', edgecolor='#2b5c8f'))

ax.set_xlabel(r'$\mathrm{Mg}^{2+}$ Reaction Coordinate Along 7.4 $\mathrm{\AA}$ Polymer Channel ($\mathrm{\AA}$)', fontweight='bold')
ax.set_ylabel(r'Relative Migration Energy (kcal/mol)', fontweight='bold')
ax.set_title(r'Cooperative Amide Bridge Lowers $\mathrm{Mg}^{2+}$ Hopping Barrier', fontweight='bold', pad=12)
ax.set_ylim(-3, 53)
ax.grid(True, linestyle=':', alpha=0.5)
ax.legend(loc='upper left', frameon=True, fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "fig2_mg_migration_barrier.png", dpi=300)
fig.savefig(OUT / "fig2_mg_migration_barrier.svg")
plt.close(fig)
print("Saved 300-DPI PNG and SVG poster figures to poster_figures/")
