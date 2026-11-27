# Nanobiohybrid Artificial Interphases for Rechargeable Magnesium-Ion Batteries (RMBs)

## Project Overview
This computational research & engineering (R&E) study investigates **biomimetic and nanobiohybrid artificial solid-electrolyte interphases (a-SEIs)** designed to lower the divalent $\text{Mg}^{2+}$ desolvation penalty in aprotic ether electrolytes (1,2-dimethoxyethane, DME). By comparing alcohol, amide, and zwitterionic coordination motifs against $\text{Mg(DME)}_3^{2+}$, we identify functional groups capable of **exergonic DME displacement** while preserving cathodic stability against Mg metal.

---

## 1. Validated DME Displacement Thermodynamics (`Figure 1`)
Replacing one bidentate DME ligand from $\text{Mg(DME)}_3^{2+}$ was evaluated across both semiempirical tight-binding (`GFN2-xTB / ALPB`) and density functional theory (`B3LYP-D3BJ / TZVP / PCM`). Both methods agree within $\sim 4\text{ kcal/mol}$:

| Candidate Ligand | Functional Motif | Resulting CN | GFN2-xTB (ALPB) | DFT B3LYP/TZVP (PCM) | Thermodynamic Outcome |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Methoxyethanol** | Monodentate Alcohol | 5 | $+23.41\text{ kcal/mol}$ | $+19.02\text{ kcal/mol}$ | Strongly endergonic (under-coordinated) |
| **DMA (`DMA_amide`)** | Monodentate Amide | 5 | $+10.77\text{ kcal/mol}$ | $+9.66\text{ kcal/mol}$ | Endergonic ($\text{CN}=6 \rightarrow 5$ penalty) |
| **Betaine** | Zwitterionic Carboxylate | **6** | **$-12.35\text{ kcal/mol}$** | **$-9.57\text{ kcal/mol}$** | **Exergonic DME stripping ($\text{CN}=6$ preserved)** |

* **Key Finding:** Zwitterionic betaine spontaneously displaces DME ($\Delta E_{\text{DFT}} = -9.57\text{ kcal/mol}$) while maintaining a six-coordinate $\text{Mg}^{2+}$ shell, avoiding the coordination-number drop ($\text{CN} = 5$) that penalizes neutral amides and alcohols.

---

## 2. Relative Vertical Reduction Ranking (`Figure 2`)
To evaluate cathodic stability, vertical electron affinities ($\text{VEA}$) were computed via unrestricted DFT (`UKS-B3LYP`) in implicit solvent and benchmarked against **ethylene carbonate (EC)** and **free DME** (`Relative Ranking, Vertical, SVP`):

| Species | Level of Theory | Vertical EA (VEA) | Est. $E_{\text{red}}$ vs. $\text{Mg}/\text{Mg}^{2+}$ |
| :--- | :--- | :---: | :---: |
| `Free_DME (Ether Solvent Ref)` | SVP (Implicit Solvent) | -1.510 eV | -3.580 V |
| `Free_DMA (Amide Motif)` | SVP (Implicit Solvent) | -0.820 eV | -2.890 V |
| `Free_Betaine (Zwitterion Motif)` | SVP (Implicit Solvent) | -0.450 eV | -2.520 V |
| `EC_Ref (Carbonate Benchmark)` | SVP (Implicit Solvent) | +1.050 eV | -1.020 V |

---

## 3. Data Hygiene & Ongoing Computational Work
* **JSON & Spin-Localization Cleanup (`results/verified_v2_data.json`):** Uncalibrated semiempirical xTB vertical electron affinities and zero-valued spin-parser artifacts are marked `not_determined` and excluded from poster figures in favor of solvent-corrected DFT benchmarks.
* **Migration Barriers (Archived Pending NEB/IRC):** Preliminary rigid forward and reverse $\text{Mg}^{2+}$ hopping scans exhibited path-dependence hysteresis and have been moved to `poster_figures/archive_pending_neb/` until full Nudged Elastic Band (NEB) and Intrinsic Reaction Coordinate (IRC) minimum-energy pathways are converged.

## Repository Contents
* `poster_figures/fig1_dme_displacement.png` (`.svg`) — GFN2-xTB vs. DFT/TZVP/PCM DME displacement energies and coordination numbers.
* `poster_figures/fig2_dft_reduction_ranking.png` (`.svg`) — Relative vertical reduction ranking benchmarked against EC and DME.
* `results/verified_v2_data.json` — Sanitized master dataset backing all poster tables and figures.
