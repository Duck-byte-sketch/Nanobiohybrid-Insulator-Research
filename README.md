# Nanobiohybrid Artificial Interphases for Rechargeable Magnesium-Ion Batteries (RMBs)

## Project Overview
This computational research & engineering (R&E) project investigates **biomimetic and nanobiohybrid artificial solid-electrolyte interphases (a-SEIs)** designed to address the divalent $\text{Mg}^{2+}$ desolvation bottleneck in aprotic ether electrolytes (1,2-dimethoxyethane, DME). By comparing alcohol, amide, and zwitterionic coordination motifs against $\text{Mg(DME)}_3^{2+}$, we evaluate two core design criteria:
1. **DME Displacement Energetics ($\Delta E$):** Whether biomimetic functional groups can exothermically strip coordinated DME without under-coordinating $\text{Mg}^{2+}$.
2. **Relative Cathodic Reducibility:** How candidate interphase ligands rank against ethylene carbonate (EC) and DME under vertical electron attachment.

---

## 1. DME Displacement Energetics (`Figure 1`)
Replacing one bidentate DME ligand from $\text{Mg(DME)}_3^{2+}$ was evaluated via electronic energies ($\Delta E$, without thermal or entropic corrections) across semiempirical tight-binding (`GFN2-xTB / ALPB`) and density functional theory (`B3LYP-D3BJ / TZVP / PCM`). Both methods agree within $\sim 4\text{ kcal/mol}$:

| Candidate Ligand | Functional Motif | Resulting CN | GFN2-xTB (ALPB) | DFT B3LYP/TZVP (PCM) | Energetic Outcome ($\Delta E$) |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Methoxyethanol** | Monodentate Alcohol | 5 | $+23.41\text{ kcal/mol}$ | $+19.02\text{ kcal/mol}$ | Strongly endothermic (under-coordinated) |
| **DMA (`DMA_amide`)** | Monodentate Amide | 5 | $+10.77\text{ kcal/mol}$ | $+9.66\text{ kcal/mol}$ | Endothermic ($\text{CN}=6 \rightarrow 5$ penalty) |
| **Betaine** | Zwitterionic Carboxylate | **6** | **$-12.35\text{ kcal/mol}$** | **$-9.57\text{ kcal/mol}$** | **Exothermic DME stripping ($\text{CN}=6$ preserved)** |

* **Key Finding:** Zwitterionic betaine exothermically displaces DME ($\Delta E_{\text{DFT}} = -9.57\text{ kcal/mol}$) while maintaining a six-coordinate $\text{Mg}^{2+}$ solvation shell, avoiding the coordination-number drop ($\text{CN} = 5$) that penalizes neutral amides and alcohols.

---

## 2. Relative Vertical Reducibility Ranking (`Figure 2`)
Vertical electron affinities ($\text{VEA}$) were computed at the **`RKS/UKS B3LYP-D3BJ / ma-def2-TZVP / IEF-PCM(THF, eps=7.58)`** level and benchmarked against **ethylene carbonate (EC)** and **free DME**:

| Species | Vertical EA (VEA) | Rel. Potential vs. $\text{Mg}/\text{Mg}^{2+}$ |
| :--- | :---: | :---: |
| `EC_Ref_Calib` | +0.576 eV | -1.49 V |
| `Free_Betaine_Calib` | +0.571 eV | -1.50 V |
| `Free_DMA_Calib` | +0.477 eV | -1.59 V |
| `Free_DME_Calib` | +0.247 eV | -1.82 V |

### Interpretation & Methodological Caveats
* **Betaine vs. EC Reducibility:** Free betaine ($+0.571\text{ eV}$) lies within $0.005\text{ eV}$ of ethylene carbonate ($+0.576\text{ eV}$)—well inside the method's uncertainty—with DMA ($+0.477\text{ eV}$) close behind and DME ($+0.247\text{ eV}$) hardest to reduce. Because EC is known to reductively passivate Mg metal, betaine is **not distinguishable from EC** under this metric, and **cathodic stability against Mg metal is not established**.
* **Relative Ranking Only:** Vertical electron affinities without adiabatic geometry relaxation or explicit cation coordination provide a **relative ordering** across ligands rather than absolute thermodynamic reduction potentials.
* **Diffuse-Basis + PCM Cavity Effect:** The small positive VEA for free DME ($+0.247\text{ eV}$, compared to $-1.51\text{ eV}$ at `def2-SVP`) is likely inflated by stabilization of diffuse excess-electron density within the continuum PCM cavity rather than representing a bound valence anion.
* **Spin Localization:** Mulliken atomic spin populations at `ma-def2-TZVP` have been scrubbed from both `results/uks_reduction_benchmark.json` and `results/verified_v2_data.json` due to basis-set ill-conditioning with diffuse functions; Löwdin/IAO populations or real-space spin densities are planned for future analysis.

---

## 3. Preliminary Constrained Hopping Scans & Next Steps
* **Constrained Multistart Barrier Scans (Archived):** Preliminary `GFN2-xTB` multistart scans with frozen $\text{Mg}$ coordinates suggest that introducing an amide bridge lowers the constrained migration envelope from $\sim 45.9\text{ kcal/mol}$ (pure betaine) to $\sim 25.3\text{ kcal/mol}$ ($\sim 45\%$ relative reduction) without creating a deep thermodynamic trap at the endpoints. However, these figures are withheld from the poster (`poster_figures/archive_pending_neb/`) because:
  1. `xtb` path optimization returned non-zero exit codes (no converged NEB/IRC saddle point yet).
  2. Multi-seed spread ($\sim 80\text{ kcal/mol}$ across 4 seeds/point) indicates upper-bound energetics rather than confirmed global minima.
  3. An absolute barrier of $\sim 25\text{ kcal/mol}$ ($\sim 1.1\text{ eV}$) remains too high for room-temperature $\text{Mg}^{2+}$ conduction without flexible backbone relaxation and DFT single-point validation.

## Repository Contents
* `poster_figures/fig1_dme_displacement.png` (`.svg`) — GFN2-xTB vs. DFT/TZVP/PCM DME displacement energies and coordination numbers.
* `poster_figures/fig2_dft_reduction_ranking.png` (`.svg`) — Relative vertical reducibility ranking at `RKS/UKS B3LYP-D3BJ / ma-def2-TZVP / IEF-PCM(THF, eps=7.58)`.
* `results/verified_v2_data.json` & `results/uks_reduction_benchmark.json` — Verified JSON datasets backing all poster tables and figures.
