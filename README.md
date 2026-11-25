# Nanobiohybrid Insulator Design for Magnesium-Ion Batteries
**Computational R&E Project (GFN2-xTB / ALPB-THF & Direct-SCF B3LYP-D3BJ/def2-SVP)**

## Overview
Multivalent $\text{Mg}^{2+}$ batteries suffer from sluggish desolvation of tightly chelated ether electrolytes (such as $\text{[Mg(DME)}_3\text{]}^{2+}$) and parasitic passivation at the Mg-metal anode. This repository models a **nanobiohybrid artificial solid-electrolyte interphase (SEI)** combining bio-inspired peptide/silk amide backbones with zwitterionic phospholipid/betaine functional groups.

## Key Results Summary

### Phase 1: Solvation Shell & Method Calibration (`AUDIT_RECORD.txt`, `pyscf_validate.py`)
* **Coordination Numbers:** $\text{[Mg(DME)}_n\text{]}^{2+}$ ($n=4,5,6$) maintains an octahedral inner shell ($\text{CN}=6$, $\text{Mg–O} = 1.93\text{–}2.03\text{ \AA}$) with outer-shell DMEs at $3.77\text{–}4.19\text{ \AA}$. Monodentate $\text{[Mg(Urea)}_n\text{]}^{2+}$ adopts $\text{CN}=4,5,6$.
* **DFT Validation (`B3LYP-D3BJ/def2-SVP`):** Total gas-phase formation energies are **$-352.57\text{ kcal/mol}$** for $\text{[Mg(DME)}_4\text{]}^{2+}$ and **$-196.70\text{ kcal/mol}$** for $\text{[Mg(Urea)}_4\text{]}^{2+}$. GFN2-xTB reproduces the qualitative ligand-exchange preference ($+84.33$ vs. $+155.88\text{ kcal/mol}$).

### Phase 2 & 3: DME-Displacement Screening (`BIOHYBRID_SCREEN.txt`, `screen_biohybrid.py`, `pyscf_verify_betaine.py`)
Reaction: $[\text{Mg(DME)}_3]^{2+} + \text{BioMotif} \longrightarrow [\text{Mg(DME)}_2(\text{BioMotif})]^{2+/+} + \text{DME}$

| Functional Motif | Bio-Polymer Analogue | CN | $\Delta E_{\text{swap}}$ (GFN2-xTB / ALPB-THF) | $\Delta E_{\text{swap}}$ (B3LYP-D3BJ/def2-SVP) | Role |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **DME (Reference)** | Baseline Ether Solvent | 6 | $0.00\text{ kcal/mol}$ | — | Symmetric Tris-Chelate |
| **Alkoxy-OH** | Cellulose / Chitosan | 5 | $+23.41\text{ kcal/mol}$ | — | Too weak to strip DME |
| **Amide (`NMA`)** | Silk Fibroin / Peptide | 5 | $+11.69\text{ kcal/mol}$ | — | Endothermic alone; acts as bulk bridge |
| **Zwitterion (`Betaine`)** | Phospholipid / Betaine | 6 | **$-12.34\text{ kcal/mol}$** | **$-36.92\text{ kcal/mol}$** | **Optimal interfacial DME stripper** |
| **Carboxylate (`OAc⁻`)** | Alginate / Pectin | 6 | $-31.14\text{ kcal/mol}$ | — | Strong desolvation; risk of trapping |

### Phase 4: Electronic Insulation & Cooperative Ion Hopping (`PHASE4_RESULTS.txt`, `phase4_transport_insulation.py`, `verify_wide_channel.py`)
* **Electronic Insulation (`B3LYP-D3BJ/def2-SVP`):** $\text{[Mg(DME)}_2\text{(Bet)]}^{2+}$ exhibits a wide **$7.82\text{ eV}$** HOMO–LUMO band gap ($\text{HOMO} = -14.10\text{ eV}$, $\text{LUMO} = -6.28\text{ eV}$), suppressing anode electron leakage.
* **Cooperative Bridge Transport (`GFN2-xTB / ALPB-THF`):** Across a $7.4\text{ \AA}$ polymer spacing between zwitterionic anchors, inserting an intermediate **amide $\text{C=O}$ bridge** reduces the $\text{Mg}^{2+}$ migration barrier by **$-20.95\text{ kcal/mol}$** ($46.38 \rightarrow 25.43\text{ kcal/mol}$), while tightly spaced ($5.2\text{ \AA}$) zwitterionic pockets achieve a **$1.76\text{ kcal/mol}$** intrinsic barrier.

## Repository Structure
* `audit.py`, `build_and_run.py`, `fix_dme_and_parse.py`, `AUDIT_RECORD.txt`: Phase 1 GFN2-xTB coordination and stepwise desolvation benchmarks.
* `pyscf_validate.py`: Memory-safe Direct-SCF `B3LYP-D3BJ/def2-SVP` benchmark script.
* `screen_biohybrid.py`, `pyscf_verify_betaine.py`, `BIOHYBRID_SCREEN.txt`: Phase 2 & 3 biohybrid ligand-displacement screen and DFT verification.
* `phase4_transport_insulation.py`, `verify_wide_channel.py`, `PHASE4_RESULTS.txt`: Phase 4 HOMO–LUMO electronic gap and two-site $\text{Mg}^{2+}$ migration barrier scans.
